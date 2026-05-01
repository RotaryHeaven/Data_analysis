%% ========================================================================
% LFP trial segmentation using NIDQ markers + Expo + TPrime
%
% Workflow:
% 1) Read NIDQ markers in NIDQ timebase
% 2) Read Expo XML and define trial structure
% 3) For each probe:
%       a) Ensure LF exists; if not, run CatGT to generate LF
%       b) Use TPrime to map NIDQ marker times to this probe timeline
%       c) Convert mapped marker times to LF sample indices
%       d) Read continuous LF trial windows directly from lf.bin
%       e) Convert raw int16 LFP to uV
%       f) Save trial data + metadata
%
% Output per probe:
%   lf_trial_data.mat
%
% Saved variables include:
%   Lfp_output_trial{block}.stim_tag
%   Lfp_output_trial{block}.stim_duration_s      -> [nTrial x 1]
%   Lfp_output_trial{block}.trial_data{trial}    -> [nKeepChan x nSamples] int16
%   Lfp_output_trial{block}.trial_data_uV{trial} -> [nKeepChan x nSamples] single
%   Lfp_output_trial{block}.trial_time{trial}    -> [1 x nSamples], stim_start = 0
%   Lfp_output_trial{block}.probe_map.acq_channel_0based
%   Lfp_output_trial{block}.probe_map.depth_um
% ========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

prestim_t  = 0.10;
poststim_t = 0.10;

root_folder = 'I:\np_data';
runName     = 'RafiL001p0120';
runind      = 1;
probes      = [0,1];     % e.g. [0 1]
expo_folder = 'I:\expo_data';

stimTag = { ...
    '[RFG_coarse2dg_99_4_150isi]', ...
    '[dir12_gpl_2_200isi_fixedphase]', ...
    '_2[Gpl2_2c_2sz_400_2_200isi]'};

% Reminder only; do not remove these channels here
needignore_site = cell(1, max(probes) + 1);
needignore_site{0 + 1} = [191];
if max(probes) >= 1
    needignore_site{1 + 1} = [14,30,191];
end

catgtExe = 'D:\workfolder_adamlab\02_Projects\awake_codes\CatGT-win\CatGT.exe';
tprimeExe = 'D:\workfolder_adamlab\02_Projects\awake_codes\TPrime-win\TPrime.exe';
sync_period_sec = 1.0;

force_rerun_event_extract = false;
force_rerun_catgt_lf      = false;
force_rerun_tprime        = false;
force_overwrite_output    = true;

manual_nidq_sync_txt = '';

maxProbeIdx = max(probes);
manual_ap_sync_txt = cell(1, maxProbeIdx + 1);
manual_ap_sync_txt(:) = {''};

%% ----------------------- Build shared session paths -----------------------

q            = @(s) ['"' s '"'];
run_g        = sprintf('%s_g%d', runName, runind);
destRoot     = fullfile(root_folder, run_g);
catgt_folder = fullfile(destRoot, ['catgt_' run_g]);
nidq_file    = fullfile(destRoot, [run_g '_t0.nidq.bin']);

if ~exist(destRoot, 'dir')
    mkdir(destRoot);
end

if ~exist(catgt_folder, 'dir')
    mkdir(catgt_folder);
end

expo_file = cell(1, numel(stimTag));
for i = 1:numel(stimTag)
    expo_file{i} = fullfile(expo_folder, sprintf('%s%s.xml', run_g, stimTag{i}));
end

fprintf('destRoot     : %s\n', destRoot);
fprintf('catgt_folder : %s\n', catgt_folder);
fprintf('nidq_file    : %s\n', nidq_file);

%% ----------------------- Basic checks -----------------------

if exist('ReadExpoXML', 'file') ~= 2
    error('ReadExpoXML.m not found on MATLAB path.');
end

if ~isfile(nidq_file)
    error('nidq_file does not exist: %s', nidq_file);
end

if ~isfile(catgtExe)
    error('CatGT executable does not exist: %s', catgtExe);
end

if ~isfile(tprimeExe)
    error('TPrime executable does not exist: %s', tprimeExe);
end

%% ----------------------- Read niSampRate from .nidq.meta -----------------------

nidq_meta_file = find_nidq_meta_file(catgt_folder, destRoot, run_g);
sr_nidq = get_meta_numeric_value(nidq_meta_file, 'niSampRate');

fprintf('nidq_meta_file: %s\n', nidq_meta_file);
fprintf('niSampRate    : %.10g\n', sr_nidq);

%% ----------------------- Extract NIDQ markers once -----------------------

done_eventunit_time = fullfile(catgt_folder, 'eventunit_time.mat');

if isfile(done_eventunit_time) && ~force_rerun_event_extract
    fprintf('\nLoading existing eventunit_time:\n  %s\n', done_eventunit_time);
    load(done_eventunit_time, 'eventunit_time');
else
    fprintf('\nExtracting NIDQ events from bin...\n');
    eventunit_time = get_eventst_GLX(nidq_file, sr_nidq);
    save(done_eventunit_time, 'eventunit_time');
    fprintf('Saved:\n  %s\n', done_eventunit_time);
end

marker_times_nidq_s = double(eventunit_time(:, 2));
nAllMarkers = numel(marker_times_nidq_s);

%% ----------------------- Read Expo and get trial definitions -----------------------

stimindex_persexpo = cell(1, numel(expo_file));
event_perexpo      = zeros(1, numel(expo_file));

for i = 1:numel(expo_file)

    if ~isfile(expo_file{i})
        error('Expo XML file does not exist: %s', expo_file{i});
    end

    expo_data = ReadExpoXML(expo_file{i});

    [stimindex_persexpo{i}, eventtime_this, ~, ~] = ...
        getExpo_pulses_fix(expo_data, 'Fixation on', 'Reward', ...
        'Ns fixation', 'fixation-start', 1);

    event_perexpo(i) = size(eventtime_this, 2);

    fprintf('Expo block %d: %s\n', i, expo_file{i});
    fprintf('  #stim trials = %d\n', numel(stimindex_persexpo{i}));
    fprintf('  #sync events = %d\n', event_perexpo(i));
end

%% ----------------------- Compute block boundaries -----------------------

cumCounts = [0, cumsum(event_perexpo)];
shift = nAllMarkers - cumCounts(end);

if cumCounts(end) > nAllMarkers
    warning('Expo sync count exceeds detected NIDQ marker count.');
end

occStart = zeros(1, numel(event_perexpo));
occEnd   = zeros(1, numel(event_perexpo));

for k = 1:numel(event_perexpo)
    if k == 3
        occStart(k) = cumCounts(k)   + 1 + shift;
        occEnd(k)   = cumCounts(k+1) + shift;
    else
        occStart(k) = cumCounts(k)   + 1;
        occEnd(k)   = cumCounts(k+1);
    end

    if occStart(k) < 1 || occEnd(k) > nAllMarkers || occStart(k) > occEnd(k)
        error('Invalid occurrence range for block %d: [%d, %d] of %d markers.', ...
            k, occStart(k), occEnd(k), nAllMarkers);
    end
end

%% ----------------------- Process each probe -----------------------

for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr   = sprintf('imec%d', thisProbe);
    this_needignore_site = needignore_site{thisProbe + 1};

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('============================================================\n');

    %% ---- Part 1: ensure LF exists ----

    lf_bin_file = fullfile(destRoot, ...
        ['catgt_' run_g], ...
        [run_g '_' imecStr], ...
        [run_g '_tcat.' imecStr '.lf.bin']);

    lf_meta_file = [lf_bin_file(1:end-3) 'meta'];
    probe_folder = fullfile(destRoot, ['catgt_' run_g], [run_g '_' imecStr]);

    if ~isfile(lf_bin_file) || force_rerun_catgt_lf

        fprintf('\nLF file not found or rerun requested. Running CatGT for probe %d...\n', thisProbe);

        params = {
            ['-dir='  q(root_folder)]
            ['-run='  runName]
            ['-g='    num2str(runind)]
            '-t=0,0'
            ['-prb='  num2str(thisProbe)]
            '-prb_fld'
            '-lf'
            '-lffilter=butter,12,0,300'
            ['-dest=' q(catgt_folder)]
            '-no_catgt_fld'
            '-out_prb_fld'
        };

        [status, cmdout, logFile] = run_catgt(catgtExe, params);
        disp(logFile)

        if status ~= 0
            error('CatGT failed for probe %d.\nOutput:\n%s', thisProbe, cmdout);
        end
    else
        fprintf('Existing LF file found:\n  %s\n', lf_bin_file);
    end

    if ~isfolder(probe_folder)
        error('probe_folder does not exist after CatGT: %s', probe_folder);
    end

    if ~isfile(lf_bin_file)
        error('LF bin still missing after CatGT: %s', lf_bin_file);
    end

    if ~isfile(lf_meta_file)
        error('LF meta still missing after CatGT: %s', lf_meta_file);
    end

    %% ---- Part 2: read LF metadata and build scale/map ----

    sr_lf       = get_meta_numeric_value(lf_meta_file, 'imSampRate');
    nSavedChans = round(get_meta_numeric_value(lf_meta_file, 'nSavedChans'));
    nKeepChans  = min(384, nSavedChans);
    nTimeLF     = get_n_timepoints_from_bin(lf_bin_file, nSavedChans);

    fprintf('lf_bin_file  : %s\n', lf_bin_file);
    fprintf('lf_meta_file : %s\n', lf_meta_file);
    fprintf('LF samp rate : %.10g Hz\n', sr_lf);
    fprintf('nSavedChans  : %d\n', nSavedChans);
    fprintf('nKeepChans   : %d\n', nKeepChans);
    fprintf('nTimeLF      : %d samples\n', nTimeLF);

    ap_bin_file  = fullfile(destRoot, ['catgt_' run_g], [run_g '_' imecStr], ...
        [run_g '_tcat.' imecStr '.ap.bin']);
    ap_meta_file = [ap_bin_file(1:end-3) 'meta'];

    meta_for_probe = choose_probe_meta_file(ap_meta_file, lf_meta_file);
    [uV_per_count, probe_map] = get_probe_scale_and_map(meta_for_probe, nKeepChans);

    %% ---- Part 3: find sync txt files for TPrime ----

    this_manual_ap_sync = manual_ap_sync_txt{thisProbe + 1};
    ap_sync_txt   = find_probe_ap_sync_txt(probe_folder, imecStr, this_manual_ap_sync);
    nidq_sync_txt = find_nidq_sync_txt(catgt_folder, manual_nidq_sync_txt);

    fprintf('AP sync txt  : %s\n', ap_sync_txt);
    fprintf('NIDQ sync txt: %s\n', nidq_sync_txt);

    %% ---- Part 4: run TPrime ----

    tprime_in_txt  = fullfile(probe_folder, 'tprime_input_nidq_marker_times.txt');
    tprime_out_txt = fullfile(probe_folder, 'tprime_output_marker_times_on_probe.txt');

    write_times_txt(tprime_in_txt, marker_times_nidq_s);

    if ~isfile(tprime_out_txt) || force_rerun_tprime
        [status, cmdout] = run_tprime( ...
            tprimeExe, sync_period_sec, ...
            ap_sync_txt, nidq_sync_txt, ...
            tprime_in_txt, tprime_out_txt);

        if status ~= 0
            error('TPrime failed for probe %d.\nOutput:\n%s', thisProbe, cmdout);
        end
    else
        fprintf('Using existing TPrime output:\n  %s\n', tprime_out_txt);
    end

    marker_times_probe_s = read_times_txt(tprime_out_txt);

    if numel(marker_times_probe_s) ~= nAllMarkers
        error(['TPrime output marker count mismatch for probe %d.\n' ...
               'Expected %d, got %d.\nFile: %s'], ...
               thisProbe, nAllMarkers, numel(marker_times_probe_s), tprime_out_txt);
    end

    %% ---- Part 5: segment LF into trials ----

    Lfp_output_trial = cell(numel(stimTag), 1);

    for j = 1:numel(stimTag)

        block_event_times_probe = marker_times_probe_s(occStart(j):occEnd(j));
        this_stim_idx = stimindex_persexpo{j};
        nTrial = numel(this_stim_idx);

        Lfp_output_trial{j} = struct();

        Lfp_output_trial{j}.stim_tag        = stimTag{j};
        Lfp_output_trial{j}.nTrial          = nTrial;
        Lfp_output_trial{j}.prestim_t       = prestim_t;
        Lfp_output_trial{j}.poststim_t      = poststim_t;
        Lfp_output_trial{j}.probe_idx       = thisProbe;
        Lfp_output_trial{j}.sr_lf           = sr_lf;
        Lfp_output_trial{j}.sr_nidq         = sr_nidq;
        Lfp_output_trial{j}.nKeepChans      = nKeepChans;
        Lfp_output_trial{j}.needignore_site = this_needignore_site(:)';
        Lfp_output_trial{j}.probe_map       = probe_map;

        Lfp_output_trial{j}.stim_duration_s = nan(nTrial, 1);
        Lfp_output_trial{j}.trial_data      = cell(nTrial, 1);
        Lfp_output_trial{j}.trial_data_uV   = cell(nTrial, 1);
        Lfp_output_trial{j}.trial_time      = cell(nTrial, 1);

        fprintf('\nProbe %d, block %d, stimTag = %s\n', thisProbe, j, stimTag{j});
        fprintf('  #trials = %d\n', nTrial);

        for m = 1:nTrial

            stim_idx = this_stim_idx(m);

            if stim_idx < 1 || (stim_idx + 1) > numel(block_event_times_probe)
                error('Invalid stim index at block %d trial %d: stim_idx = %d', ...
                    j, m, stim_idx);
            end

            stim_start_probe = block_event_times_probe(stim_idx);
            stim_end_probe   = block_event_times_probe(stim_idx + 1);
            stim_duration_s  = stim_end_probe - stim_start_probe;

            trial_t0_probe = stim_start_probe - prestim_t;
            trial_t1_probe = stim_end_probe   + poststim_t;

            req_sample0 = floor(trial_t0_probe * sr_lf) + 1;
            req_sample1 = ceil(trial_t1_probe * sr_lf)  + 1;

            [trial_data, sample_idx_read, ~, ~, ~] = ...
                read_lf_window(lf_bin_file, nSavedChans, nKeepChans, nTimeLF, req_sample0, req_sample1);

            trial_time_rel = ((double(sample_idx_read) - 1) ./ sr_lf) - stim_start_probe;

            Lfp_output_trial{j}.stim_duration_s(m)  = stim_duration_s;
            Lfp_output_trial{j}.trial_data{m}       = trial_data;
            Lfp_output_trial{j}.trial_data_uV{m}    = single(trial_data) .* single(uV_per_count(:));
            Lfp_output_trial{j}.trial_time{m}       = trial_time_rel(:)';
        end
    end

    %% ---- Part 6: save output for this probe ----

    out_file = fullfile(probe_folder, 'lf_trial_data.mat');

    if isfile(out_file) && ~force_overwrite_output
        warning('Output exists and overwrite disabled, skipping save: %s', out_file);
    else
        save(out_file, ...
            'Lfp_output_trial', ...
            'stimTag', ...
            'prestim_t', ...
            'poststim_t', ...
            'runName', ...
            'runind', ...
            'run_g', ...
            '-v7.3');

        fprintf('\nSaved:\n  %s\n', out_file);
    end
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================

function [status, cmdout, logFile] = run_catgt(catgtExe, paramList)
assert(isfile(catgtExe), 'CatGT.exe does not exist: %s', catgtExe);

dest = '';
for k = 1:numel(paramList)
    tok = regexp(paramList{k}, '^-dest=(.+)$', 'tokens', 'once');
    if ~isempty(tok)
        dest = tok{1};
        break;
    end
end
assert(~isempty(dest), 'paramList needs -dest=...');

args = strjoin(paramList, ' ');
cmd  = sprintf('pushd "%s" && "%s" %s && popd', dest, catgtExe, args);

[status, cmdout] = system(cmd);
logFile = fullfile(dest, 'catGT.log');

if status == 0
    fprintf('[CatGT] success\n');
else
    warning('[CatGT] failed, status=%d.\nOutput:\n%s', status, cmdout);
end
end


function meta_file = find_nidq_meta_file(catgt_folder, destDir, run_g)
meta_list = dir(fullfile(catgt_folder, '*.nidq.meta'));

if isempty(meta_list)
    meta_list = dir(fullfile(destDir, '*.nidq.meta'));
end

if isempty(meta_list)
    error('No .nidq.meta file found in:\n%s\nor\n%s', catgt_folder, destDir);
end

if numel(meta_list) == 1
    meta_file = fullfile(meta_list(1).folder, meta_list(1).name);
    return;
end

fprintf('\nMultiple .nidq.meta files were found for run %s:\n', run_g);
for i = 1:numel(meta_list)
    fprintf('  [%d] %s\n', i, fullfile(meta_list(i).folder, meta_list(i).name));
end

while true
    idx = input(sprintf('Please choose one meta file [1-%d]: ', numel(meta_list)));
    if isempty(idx)
        fprintf('No input detected. Please enter a number.\n');
        continue;
    end

    if ~isscalar(idx) || ~isnumeric(idx) || isnan(idx) || ...
            idx < 1 || idx > numel(meta_list) || floor(idx) ~= idx
        fprintf('Invalid choice. Please enter an integer from 1 to %d.\n', numel(meta_list));
        continue;
    end
    break;
end

meta_file = fullfile(meta_list(idx).folder, meta_list(idx).name);
fprintf('Selected meta file: %s\n', meta_file);
end


function val = get_meta_numeric_value(meta_file, key)
if ~isfile(meta_file)
    error('Meta file does not exist: %s', meta_file);
end

fid = fopen(meta_file, 'r');
if fid == -1
    error('Cannot open meta file: %s', meta_file);
end
cleanupObj = onCleanup(@() fclose(fid));

val = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ~ischar(line) || isempty(line)
        continue;
    end

    prefix = [key '='];
    if startsWith(line, prefix)
        val_str = strtrim(line(numel(prefix)+1:end));
        val = str2double(val_str);
        if isnan(val)
            error('Meta value for %s is not numeric in %s: %s', key, meta_file, val_str);
        end
        return;
    end
end

error('Key %s not found in meta file: %s', key, meta_file);
end


function eventunit_time = get_eventst_GLX(fn, sr_nidq)
if exist(fn, 'file') ~= 2
    [filename, loadpath] = uigetfile('*.bin', 'Please select NIDQ .bin file');
    if isequal(filename, 0)
        error('No nidq file selected.');
    end
    fn = fullfile(loadpath, filename);
end

n_chan = 4;
d_size = round(sr_nidq * 30);

fni = dir(fn);
n_sample = floor(fni(1).bytes / 2);
n_time   = floor(n_sample / n_chan);

sync_channel = zeros(1, n_time, 'int16');
t = 1;

fid = fopen(fn, 'r');
if fid == -1
    error('Cannot open file: %s', fn);
end
cleanupObj = onCleanup(@() fclose(fid));

counter = 0;
fseek(fid, 0, 'bof');

while ~feof(fid)
    d_size_t = min(d_size, floor((n_sample - counter) / n_chan));

    if d_size_t > 0
        [data, add_count] = fread(fid, [n_chan d_size_t], '*int16');
        sync_channel(1, t:t+d_size_t-1) = bitget(data(end, :), 2);
        counter = counter + add_count;
        t = t + d_size_t;
    else
        break;
    end
end

assert(counter == n_sample, 'Could not read all samples');

eventsample = find(diff(sync_channel) > 0.5);
eventt = double(eventsample(:)) / sr_nidq;

fprintf('Assumed nidq sampling rate: %.10g\n', sr_nidq);
fprintf('Detected %d events\n', numel(eventt));

eventunit_time = [ones(numel(eventt), 1) * 2000, eventt];
end


function txt_file = find_probe_ap_sync_txt(probe_folder, imecStr, manual_txt)
if ~isempty(manual_txt)
    if ~isfile(manual_txt)
        error('Manual AP sync txt does not exist: %s', manual_txt);
    end
    txt_file = manual_txt;
    return;
end

pat = fullfile(probe_folder, ['*.' imecStr '.ap.xd_*.txt']);
d = dir(pat);

if isempty(d)
    error(['No AP sync-edge txt found with pattern:' newline '%s' newline ...
           'Set manual_ap_sync_txt for this probe if needed.'], pat);
end

if numel(d) == 1
    txt_file = fullfile(d(1).folder, d(1).name);
    return;
end

fprintf('\nMultiple AP xd txt files found in %s\n', probe_folder);
for i = 1:numel(d)
    fprintf('  [%d] %s\n', i, fullfile(d(i).folder, d(i).name));
end
error('Multiple AP xd txt files found. Please set manual_ap_sync_txt for this probe.');
end


function txt_file = find_nidq_sync_txt(catgt_folder, manual_txt)
if ~isempty(manual_txt)
    if ~isfile(manual_txt)
        error('Manual NIDQ sync txt does not exist: %s', manual_txt);
    end
    txt_file = manual_txt;
    return;
end

pat = fullfile(catgt_folder, '*.nidq.xd_*.txt');
d = dir(pat);

if isempty(d)
    error(['No NIDQ sync-edge txt found with pattern:' newline '%s' newline ...
           'Set manual_nidq_sync_txt if needed.'], pat);
end

if numel(d) == 1
    txt_file = fullfile(d(1).folder, d(1).name);
    return;
end

fprintf('\nMultiple NIDQ xd txt files found in %s\n', catgt_folder);
for i = 1:numel(d)
    fprintf('  [%d] %s\n', i, fullfile(d(i).folder, d(i).name));
end
error('Multiple NIDQ xd txt files found. Please set manual_nidq_sync_txt.');
end


function write_times_txt(out_file, t)
t = double(t(:));

fid = fopen(out_file, 'w');
if fid == -1
    error('Cannot open file for writing: %s', out_file);
end
cleanupObj = onCleanup(@() fclose(fid));

for i = 1:numel(t)
    fprintf(fid, '%.9f\n', t(i));
end
end


function t = read_times_txt(in_file)
if ~isfile(in_file)
    error('Times txt file does not exist: %s', in_file);
end

fid = fopen(in_file, 'r');
if fid == -1
    error('Cannot open file for reading: %s', in_file);
end
cleanupObj = onCleanup(@() fclose(fid));

t = fscanf(fid, '%f');
t = double(t(:));

if isempty(t)
    error('No times read from file: %s', in_file);
end

if any(~isfinite(t))
    error('Non-finite times found in file: %s', in_file);
end
end


function [status, cmdout] = run_tprime(tprimeExe, syncperiod, tostream_txt, fromstream_txt, events_in_txt, events_out_txt)
assert(isfile(tprimeExe), 'TPrime.exe does not exist: %s', tprimeExe);
assert(isfile(tostream_txt), 'tostream txt missing: %s', tostream_txt);
assert(isfile(fromstream_txt), 'fromstream txt missing: %s', fromstream_txt);
assert(isfile(events_in_txt), 'events_in txt missing: %s', events_in_txt);

q = @(s) ['"' s '"'];

cmd = sprintf('%s -syncperiod=%.10g -tostream=%s -fromstream=1,%s -events=1,%s,%s', ...
    q(tprimeExe), syncperiod, q(tostream_txt), q(fromstream_txt), q(events_in_txt), q(events_out_txt));

fprintf('\nRunning TPrime:\n%s\n', cmd);

[status, cmdout] = system(cmd);

if status == 0
    fprintf('[TPrime] success\n');
else
    warning('[TPrime] failed, status=%d.\nOutput:\n%s', status, cmdout);
end
end


function nTime = get_n_timepoints_from_bin(bin_file, nSavedChans)
if ~isfile(bin_file)
    error('Bin file does not exist: %s', bin_file);
end

d = dir(bin_file);
n_int16 = d.bytes / 2;

if rem(n_int16, nSavedChans) ~= 0
    error('File size is not divisible by nSavedChans.\nFile: %s\nnSavedChans=%d', ...
        bin_file, nSavedChans);
end

nTime = n_int16 / nSavedChans;
end


function [data_keep, sample_idx, s0, s1, was_clipped] = read_lf_window(bin_file, nSavedChans, nKeepChans, nTimeTotal, req_s0, req_s1)
if req_s1 < req_s0
    error('Requested sample window is invalid: [%d, %d]', req_s0, req_s1);
end

s0 = max(1, req_s0);
s1 = min(nTimeTotal, req_s1);
was_clipped = (s0 ~= req_s0) || (s1 ~= req_s1);

if s1 < s0
    error('Requested LF window falls completely outside file bounds.');
end

nSamp = s1 - s0 + 1;
sample_idx = s0:s1;

fid = fopen(bin_file, 'r');
if fid == -1
    error('Cannot open LF bin file: %s', bin_file);
end
cleanupObj = onCleanup(@() fclose(fid));

byte_offset = (s0 - 1) * nSavedChans * 2;
fseek_status = fseek(fid, byte_offset, 'bof');
if fseek_status ~= 0
    error('fseek failed for file: %s', bin_file);
end

data = fread(fid, [nSavedChans, nSamp], '*int16');

if size(data, 2) ~= nSamp
    error('Could not read requested LF window.\nRequested %d samples, got %d.', ...
        nSamp, size(data, 2));
end

keep_idx = 1:nKeepChans;
data_keep = data(keep_idx, :);
end


function meta_file = choose_probe_meta_file(ap_meta_file, lf_meta_file)
if isfile(ap_meta_file)
    meta_file = ap_meta_file;
elseif isfile(lf_meta_file)
    meta_file = lf_meta_file;
else
    error('Neither AP nor LF meta exists.\nAP: %s\nLF: %s', ap_meta_file, lf_meta_file);
end
end


function [uV_per_count, probe_map] = get_probe_scale_and_map(meta_file, nKeepChans)
meta = ReadMeta_local(meta_file, '');

saved_ch_1based = OriginalChans_local(meta);
if numel(saved_ch_1based) < nKeepChans
    error('nKeepChans exceeds number of saved channels in meta.');
end
saved_ch_1based = saved_ch_1based(1:nKeepChans);

if isfield(meta,'snsGeomMap')
    [~, ~, ~, ~, ~, depth_all, ~] = geomMapToGeom_local(meta);
elseif isfield(meta,'snsShankMap')
    [~, ~, ~, ~, ~, depth_all, ~] = shankMapToGeom_local(meta);
else
    error('Meta lacks snsGeomMap/snsShankMap: %s', meta_file);
end

depth_all = depth_all(:);
if numel(depth_all) < nKeepChans
    error('Geometry depth array shorter than nKeepChans.');
end
depth_um = depth_all(1:nKeepChans);

fI2V = Int2VoltsIM_local(meta);
[APgain, LFgain] = ChanGainsIM_local(meta);

use_lf_gain = ~isempty(LFgain) && any(LFgain > 0);

uV_per_count = nan(nKeepChans, 1);
for i = 1:nKeepChans
    k = saved_ch_1based(i);

    if use_lf_gain && k <= numel(LFgain) && LFgain(k) > 0
        gain = LFgain(k);
    elseif k <= numel(APgain) && APgain(k) > 0
        gain = APgain(k);
    elseif k <= numel(LFgain) && LFgain(k) > 0
        gain = LFgain(k);
    else
        error('Cannot determine gain for saved channel row %d (acq idx %d).', i, k-1);
    end

    uV_per_count(i) = 1e6 * fI2V / gain;
end

probe_map = struct();
probe_map.acq_channel_0based = saved_ch_1based(:) - 1;
probe_map.depth_um           = depth_um(:);
end


function meta = ReadMeta_local(metaName, path)
fid = fopen(fullfile(path, metaName), 'r');
if fid == -1
    error('Cannot open meta file: %s', fullfile(path, metaName));
end
c = onCleanup(@() fclose(fid));

C = textscan(fid, '%[^=] = %[^\r\n]');

meta = struct();

for i = 1:length(C{1})
    tag = strtrim(C{1}{i});
    if ~isempty(tag) && tag(1) == '~'
        tag = tag(2:end);
    end
    meta.(tag) = strtrim(C{2}{i});
end
end


function fI2V = Int2VoltsIM_local(meta)
if isfield(meta, 'imMaxInt')
    maxInt = str2double(meta.imMaxInt);
else
    maxInt = 512;
end

fI2V = str2double(meta.imAiRangeMax) / maxInt;
end


function chans = OriginalChans_local(meta)
if strcmp(meta.snsSaveChanSubset, 'all')
    chans = (1:str2double(meta.nSavedChans))';
else
    chans = str2num(meta.snsSaveChanSubset)'; %#ok<ST2NM>
    chans = chans + 1;
end
end


function [APgain, LFgain] = ChanGainsIM_local(meta)
np1_imro = [0,1020,1030,1200,1100,1120,1121,1122,1123,1300];

if isfield(meta, 'acqApLfSy')
    acqCountList = str2num(meta.acqApLfSy); %#ok<ST2NM>
elseif isfield(meta, 'snsApLfSy')
    acqCountList = str2num(meta.snsApLfSy); %#ok<ST2NM>
else
    error('Meta missing acqApLfSy/snsApLfSy.');
end

nAP = acqCountList(1);
nLF = acqCountList(2);

APgain = zeros(nAP, 1);
LFgain = zeros(nLF, 1);

if isfield(meta, 'imDatPrb_type')
    probeType = str2double(meta.imDatPrb_type);
else
    probeType = 0;
end

if ismember(probeType, np1_imro)
    if isfield(meta, 'typeEnabled')
        C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
            'EndOfLine', ')', 'HeaderLines', 1);
    else
        C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
            'EndOfLine', ')', 'HeaderLines', 1);
    end

    APgain = double(cell2mat(C(1)));
    LFgain = double(cell2mat(C(2)));

else
    if isfield(meta, 'imChan0apGain')
        APgain = APgain + str2double(meta.imChan0apGain);

        if nLF > 0 && isfield(meta, 'imChan0lfGain')
            LFgain = LFgain + str2double(meta.imChan0lfGain);
        end

    elseif probeType == 1110
        currList = sscanf(meta.imroTbl, '(%d,%d,%d,%d,%d');
        APgain = APgain + currList(4);
        if nLF > 0
            LFgain = LFgain + currList(5);
        end

    elseif probeType == 21 || probeType == 24
        APgain = APgain + 80;

    elseif probeType == 2013
        APgain = APgain + 100;

    else
        warning('Unknown probe type %g, fallback gain=1.', probeType);
        APgain = APgain + 1;
        if nLF > 0
            LFgain = LFgain + 1;
        end
    end
end
end


function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom_local(meta)
C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );

shankInd = double(cell2mat(C(1)));
xCoord   = double(cell2mat(C(2)));
yCoord   = double(cell2mat(C(3)));
connected = double(cell2mat(C(4)));

geomStr = meta.snsGeomMap;
headStr = extractBefore(geomStr,')(');
headParts = split(headStr,',');
nShank     = str2double(headParts{2});
shankWidth = str2double(headParts{4});
shankPitch = str2double(headParts{3});
end


function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom_local(meta)
[nchan,~,~] = ChannelCountsIM_local(meta);

C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );

shankInd = double(cell2mat(C(1)));
colInd   = double(cell2mat(C(2)));
rowInd   = double(cell2mat(C(3)));
connected = double(cell2mat(C(4)));

shankInd  = shankInd(1:nchan);
colInd    = colInd(1:nchan);
rowInd    = rowInd(1:nchan);
connected = connected(1:nchan);

geom = getGeomParams_local(meta);

oddRows  = logical(mod(rowInd,2));
evenRows = ~oddRows;

xCoord = colInd * geom.horzPitch;
xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff;
xCoord(oddRows)  = xCoord(oddRows)  + geom.odd_xOff;
yCoord = rowInd * geom.vertPitch;

nShank     = geom.nShank;
shankWidth = geom.shankWidth;
shankPitch = geom.shankPitch;
end


function [AP, LF, SY] = ChannelCountsIM_local(meta)
if isfield(meta, 'snsApLfSy')
    M = str2num(meta.snsApLfSy); %#ok<ST2NM>
elseif isfield(meta, 'acqApLfSy')
    M = str2num(meta.acqApLfSy); %#ok<ST2NM>
else
    error('Meta missing snsApLfSy/acqApLfSy.');
end

AP = M(1);
LF = M(2);
SY = M(3);
end


function geom = getGeomParams_local(meta)
geomTypeMap = makeTypeMap_local();

if isfield(meta,'imDatPrb_pn')
    pn = meta.imDatPrb_pn;
else
    pn = '3A';
end

if geomTypeMap.isKey(pn)
    geomType = geomTypeMap(pn);
else
    error('Unsupported probe part number: %s', pn);
end

switch geomType
    case 'np1_stag_70um'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 11;
        geom.horzPitch = 32; geom.vertPitch = 20;
    case 'nhp_lin_70um'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 27;
        geom.horzPitch = 32; geom.vertPitch = 20;
    case 'nhp_stag_125um_med'
        geom.nShank = 1; geom.shankWidth = 125; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 11;
        geom.horzPitch = 87; geom.vertPitch = 20;
    case 'nhp_stag_125um_long'
        geom.nShank = 1; geom.shankWidth = 125; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 11;
        geom.horzPitch = 87; geom.vertPitch = 20;
    case 'nhp_lin_125um_med'
        geom.nShank = 1; geom.shankWidth = 125; geom.shankPitch = 0;
        geom.even_xOff = 11; geom.odd_xOff = 11;
        geom.horzPitch = 103; geom.vertPitch = 20;
    case 'nhp_lin_125um_long'
        geom.nShank = 1; geom.shankWidth = 125; geom.shankPitch = 0;
        geom.even_xOff = 11; geom.odd_xOff = 11;
        geom.horzPitch = 103; geom.vertPitch = 20;
    case 'uhd_8col_1bank'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 14; geom.odd_xOff = 14;
        geom.horzPitch = 6; geom.vertPitch = 6;
    case 'uhd_8col_16bank'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 14; geom.odd_xOff = 14;
        geom.horzPitch = 6; geom.vertPitch = 6;
    case 'np2_ss'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 27;
        geom.horzPitch = 32; geom.vertPitch = 15;
    case 'np2_4s'
        geom.nShank = 4; geom.shankWidth = 70; geom.shankPitch = 250;
        geom.even_xOff = 27; geom.odd_xOff = 27;
        geom.horzPitch = 32; geom.vertPitch = 15;
    case 'NP1120'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 6.75; geom.odd_xOff = 6.75;
        geom.horzPitch = 4.5; geom.vertPitch = 4.5;
    case 'NP1121'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 6.25; geom.odd_xOff = 6.25;
        geom.horzPitch = 3; geom.vertPitch = 3;
    case 'NP1122'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 12.5; geom.odd_xOff = 12.5;
        geom.horzPitch = 3; geom.vertPitch = 3;
    case 'NP1123'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 10.25; geom.odd_xOff = 10.25;
        geom.horzPitch = 4.5; geom.vertPitch = 4.5;
    case 'NP1300'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 11; geom.odd_xOff = 11;
        geom.horzPitch = 48; geom.vertPitch = 20;
    case 'NP1200'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 27; geom.odd_xOff = 11;
        geom.horzPitch = 32; geom.vertPitch = 20;
    case 'NXT3000'
        geom.nShank = 1; geom.shankWidth = 70; geom.shankPitch = 0;
        geom.even_xOff = 53; geom.odd_xOff = 53;
        geom.horzPitch = 0; geom.vertPitch = 15;
    otherwise
        error('Unsupported geometry type: %s', geomType);
end
end


function M = makeTypeMap_local()
M = containers.Map('KeyType','char','ValueType','char');

M('3A') = 'np1_stag_70um';
M('PRB_1_4_0480_1') = 'np1_stag_70um';
M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
M('NP1010') = 'np1_stag_70um';
M('NP1011') = 'np1_stag_70um';
M('NP1012') = 'np1_stag_70um';
M('NP1013') = 'np1_stag_70um';

M('NP1015') = 'nhp_lin_70um';
M('NP1016') = 'nhp_lin_70um';
M('NP1017') = 'nhp_lin_70um';

M('NP1020') = 'nhp_stag_125um_med';
M('NP1021') = 'nhp_stag_125um_med';
M('NP1030') = 'nhp_stag_125um_long';
M('NP1031') = 'nhp_stag_125um_long';

M('NP1022') = 'nhp_lin_125um_med';
M('NP1032') = 'nhp_lin_125um_long';

M('NP1100') = 'uhd_8col_1bank';
M('NP1110') = 'uhd_8col_16bank';

M('PRB2_1_2_0640_0') = 'np2_ss';
M('PRB2_1_4_0480_1') = 'np2_ss';
M('NP2000') = 'np2_ss';
M('NP2003') = 'np2_ss';
M('NP2004') = 'np2_ss';

M('PRB2_4_2_0640_0') = 'np2_4s';
M('PRB2_4_4_0480_1') = 'np2_4s';
M('NP2010') = 'np2_4s';
M('NP2013') = 'np2_4s';
M('NP2014') = 'np2_4s';

M('NP1120') = 'NP1120';
M('NP1121') = 'NP1121';
M('NP1122') = 'NP1122';
M('NP1123') = 'NP1123';
M('NP1300') = 'NP1300';

M('NP1200') = 'NP1200';
M('NXT3000') = 'NXT3000';
end