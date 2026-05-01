%% ========================================================================
%% Notes
% This script processes one SpikeGLX/Kilosort session for a selected run and probe.
% First, set user parameters, folder paths, run name, probe index, and Expo stimulus tags.
% The script then builds all file paths and reads niSampRate from the .nidq.meta file.
% NIDQ event markers are extracted once per run and saved in catgt_folder/eventunit_time.mat.
% If eventunit_time.mat already exists, the event extraction step is skipped.
% After that, each kilosort* folder under the probe folder is processed independently.
% For each kilosort folder, spike times are combined with event markers, then split into
% occurrence-level blocks and trial-level windows around each stimulus.
% Run summary metrics and figures are generated for each kilosort folder.
% If all expected output files already exist in a kilosort folder, that folder is skipped.
% Errors in one kilosort folder do not stop the whole script; processing continues to the next folder.
% ========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

% Time window around each stimulus
prestim_t  = 0.1;
poststim_t = 0.1;

%% ----------------------- Directory settings -----------------------

root_folder = 'I:\np_data';
runName     = 'RafiL001p0120';
runind      = 1;         % run index after -g
probes      = [0,1];    % probe indices after -prb
expo_folder = 'I:\expo_data';

% Expo stimulus files to use for this session
stimTag = { ...
    '[RFG_coarse2dg_99_4_150isi]', ...
    '[dir12_gpl_2_200isi_fixedphase]', ...
    '_2[Gpl2_2c_2sz_400_2_200isi]'};

%% ----------------------- Build shared session paths -----------------------

run_g   = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);

% NIDQ binary file for this session
nidq_file = fullfile(destDir, [run_g '_t0.nidq.bin']);

% Expo XML files for each stimulus block
expo_file = cell(1, numel(stimTag));
for i = 1:numel(stimTag)
    expo_file{i} = fullfile(expo_folder, sprintf('%s%s.xml', run_g, stimTag{i}));
end

% CatGT folder shared by this run
catgt_folder = fullfile(destDir, ['catgt_' run_g]);

%% ----------------------- Read niSampRate from .nidq.meta -----------------------

% Find the .nidq.meta file, then extract niSampRate from it
nidq_meta_file = find_nidq_meta_file(catgt_folder, destDir, run_g);
sr_nidq = get_niSampRate(nidq_meta_file);

fprintf('destDir       : %s\n', destDir);
fprintf('catgt_folder  : %s\n', catgt_folder);
fprintf('nidq_file     : %s\n', nidq_file);
fprintf('nidq_meta_file: %s\n', nidq_meta_file);
fprintf('niSampRate    : %.10g\n', sr_nidq);

%% ----------------------- Basic checks -----------------------

% readNPY is needed to load Kilosort .npy files
if exist('readNPY', 'file') ~= 2
    error(['readNPY.m was not found in the MATLAB path. ', ...
        'Please add npy-matlab to your path first.']);
end

if ~isfile(nidq_file)
    error('nidq_file does not exist: %s', nidq_file);
end

%% ----------------------- Extract NIDQ events once -----------------------

% eventunit_time:
%   column 1 = unit id, fixed as 2000 for inserted event markers
%   column 2 = event time in seconds, in the NIDQ timebase

done_eventunit_time = fullfile(catgt_folder, 'eventunit_time.mat');

if isfile(done_eventunit_time)
    fprintf('\nSkipping expo extract event step (already analyzed): %s\n', catgt_folder);
    load(done_eventunit_time, 'eventunit_time');
else
    fprintf('\nExtracting NIDQ events once for this session...\n');
    eventunit_time = get_eventst_GLX(nidq_file, sr_nidq);

    fprintf('Saved shared eventunit_time in cat folder:\n');
    fprintf('  %s\n', done_eventunit_time);
    save(done_eventunit_time, 'eventunit_time');
end

%% ----------------------- Read Expo files and segment info -----------------------

stimindex_persexpo = cell(1, numel(expo_file));
eventtime_perexpo  = cell(1, numel(expo_file));
event_perexpo      = zeros(1, numel(expo_file));

for i = 1:numel(expo_file)

    if ~isfile(expo_file{i})
        error('Expo XML file does not exist: %s', expo_file{i});
    end

    % Read one Expo XML file
    expo_data = ReadExpoXML(expo_file{i});

    % Extract completed trial event indices and event times from Expo
    [stimindex_persexpo{i}, eventtime_perexpo{i}, ~, ~] = ...
        getExpo_pulses_fix(expo_data, 'Fixation on', 'Reward', ...
        'Ns fixation', 'fixation-start', 1);

    % Number of sync events in this Expo block
    event_perexpo(i) = size(eventtime_perexpo{i}, 2);
end

%% ----------------------- Compute segment boundaries in eventunit_time -----------------------

% cumCounts gives cumulative event counts across Expo blocks
cumCounts = [0, cumsum(event_perexpo)];

% If NIDQ has extra leading events relative to Expo, shift later indexing
shift = size(eventunit_time, 1) - cumCounts(end);

if cumCounts(end) > size(eventunit_time, 1)
    warning('Unequal number of syncs in expo and the nidq files');
end

% occStart / occEnd define the event index range for each block
occStart = zeros(1, numel(event_perexpo));
occEnd   = zeros(1, numel(event_perexpo));
dur      = cell(1, numel(event_perexpo));
expotime = cell(1, numel(event_perexpo));
nidqtime = cell(1, numel(event_perexpo));

for k = 1:numel(event_perexpo)

    % Custom alignment rule kept from your current logic:
    % block 3 is shifted by "shift", other blocks are not
    if k == 3
        occStart(k) = cumCounts(k)   + 1 + shift;
        occEnd(k)   = cumCounts(k+1) + shift;
    else
        occStart(k) = cumCounts(k)   + 1;
        occEnd(k)   = cumCounts(k+1);
    end

    % Time column for this event block
    t = eventunit_time(occStart(k):occEnd(k), 2);

    % Pairwise event-to-event duration within this block
    dur{k} = diff(t);

    % Compare Expo and NIDQ timing for sanity checking
    [expotime{k}, nidqtime{k}] = run_exponidq_sanity_check( ...
        dur{k}, eventtime_perexpo{k}, catgt_folder, sr_nidq);
end

close all

%% ----------------------- Process each probe folder -----------------------

for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr   = sprintf('imec%d', thisProbe);

    probe_folder = fullfile( ...
        destDir, ...
        ['catgt_' run_g], ...
        [run_g '_' imecStr] ...
        );

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    if ~isfolder(probe_folder)
        warning('probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    %% ----------------------- Find all kilosort folders for this probe -----------------------

    % Each kilosort* folder under the probe folder will be processed
    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);

    if isempty(d)
        warning('No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    % Sort for stable processing order
    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf('Found %d kilosort folder(s) under probe %d.\n', numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------

    for i = 1:numel(d)

        ksDir = fullfile(d(i).folder, d(i).name);

        % skip this kilosort folder if the analysis outputs already exist
        done_file_1 = fullfile(ksDir, 'spike_unit_time.mat');
        done_file_2 = fullfile(ksDir, 'spike_unit_time_occ.mat');
        done_file_3 = fullfile(ksDir, 'spike_unit_time_trial.mat');
        done_file_4 = fullfile(ksDir, 'run_metrics_summary.png');
        done_file_5 = fullfile(ksDir, 'run_metrics_summary.fig');
        done_file_6 = fullfile(ksDir, 'unit_run_metrics.mat');

        if isfile(done_file_1) && isfile(done_file_2) && isfile(done_file_3) && ...
                isfile(done_file_4) && isfile(done_file_5) && isfile(done_file_6)
            fprintf('\nSkipping ksDir (already analyzed): %s\n', ksDir);
            continue;
        end

        fprintf('\nProcessing probe %d, ksDir: %s\n', thisProbe, ksDir);

        try
            %% Build spike_unit_time for this kilosort folder

            % spike_unit_time:
            %   column 1 = unit id
            %   column 2 = spike time in seconds (already mapped to NIDQ timebase)
            spike_unit_time = build_spike_unit_time_from_ksdir(ksDir);
            Allunits = unique(spike_unit_time(:,1));

            % Append NIDQ events as unit 2000
            spike_unit_time = [spike_unit_time; eventunit_time];

            % Sort the combined matrix by time
            spike_unit_time = sortrows(spike_unit_time, 2);

            % Save the full combined matrix
            save(fullfile(ksDir, 'spike_unit_time.mat'), 'spike_unit_time');

            fprintf('Saved:\n');
            fprintf('  %s\n', fullfile(ksDir, 'spike_unit_time.mat'));

            %% Split spike_unit_time into larger occurrence blocks

            % Find the rows corresponding to the inserted event markers
            event_rows_in_spike = find(spike_unit_time(:, 1) == 2000);

            % Number of occurrence blocks
            nOcc = numel(occStart);

            % spike_unit_time_occ{j} contains one larger block between event indices
            spike_unit_time_occ = cell(nOcc, 1);

            for j = 1:nOcc

                % Convert event indices into row indices in the sorted combined matrix
                row_start = event_rows_in_spike(occStart(j));
                row_end   = event_rows_in_spike(occEnd(j));

                % Keep everything between the start and end events, inclusive
                spike_unit_time_occ{j} = spike_unit_time(row_start:row_end, :);
            end

            save(fullfile(ksDir, 'spike_unit_time_occ.mat'), 'spike_unit_time_occ');
            fprintf('  %s\n', fullfile(ksDir, 'spike_unit_time_occ.mat'));

            %% Split each occurrence block into trial-level windows

            % spike_unit_time_trial{j}{m} = m-th trial in occurrence block j
            spike_unit_time_trial = cell(numel(spike_unit_time_occ), 1);

            for j = 1:numel(spike_unit_time_occ)

                this_occ = spike_unit_time_occ{j};

                % Event rows and times within this occurrence block
                event_times = this_occ(this_occ(:, 1) == 2000, 2);

                % Expo-derived stimulus start indices for this block
                this_stim_idx = stimindex_persexpo{j};
                nTrial = numel(this_stim_idx);

                spike_unit_time_trial{j} = cell(nTrial, 1);

                for m = 1:nTrial
                    stim_idx = this_stim_idx(m);

                    % Trial start = event stim_idx
                    % Trial end   = event stim_idx + 1
                    if stim_idx < 1 || (stim_idx + 1) > numel(event_times)
                        error('Invalid stim index at occ %d trial %d: stim_idx = %d', ...
                            j, m, stim_idx);
                    end

                    stim_start = event_times(stim_idx);
                    stim_end   = event_times(stim_idx + 1);

                    % Window around the stimulus period
                    trial_t0 = stim_start - prestim_t;
                    trial_t1 = stim_end   + poststim_t;

                    % Keep all spikes/events falling inside this trial window
                    keep_idx = this_occ(:, 2) >= trial_t0 & this_occ(:, 2) <= trial_t1;

                    % Shift time so that stim_start becomes time zero
                    tmp = this_occ(keep_idx, :);
                    tmp(:, 2) = tmp(:, 2) - stim_start;

                    spike_unit_time_trial{j}{m} = tmp;
                end
            end

            save(fullfile(ksDir, 'spike_unit_time_trial.mat'), ...
                'spike_unit_time_trial', 'prestim_t', 'poststim_t');
            fprintf('  %s\n', fullfile(ksDir, 'spike_unit_time_trial.mat'));

            % plot per-run summary metrics for this kilosort folder
            unit_run_metrics = plot_kilosort_run_metrics( ...
                spike_unit_time_trial, Allunits, stimTag, ...
                prestim_t, poststim_t, ksDir);

            save(fullfile(ksDir, 'unit_run_metrics.mat'), 'unit_run_metrics');

            fprintf('  %s\n', fullfile(ksDir, 'run_metrics_summary.png'));
            fprintf('  %s\n', fullfile(ksDir, 'run_metrics_summary.fig'));
            fprintf('  %s\n', fullfile(ksDir, 'unit_run_metrics.mat'));

        catch ME
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', ME.message);
        end
    end
end

fprintf('\nDone.\n');

%% ----------------------- Local functions -----------------------

function spike_unit_time = build_spike_unit_time_from_ksdir(ksDir)
% Read one Kilosort folder and construct spike_unit_time.
%
% Output:
%   spike_unit_time(:,1) = unit id from spike_clusters.npy
%   spike_unit_time(:,2) = spike time in seconds from spike_seconds_adj.npy

if ~isfolder(ksDir)
    error('ksDir does not exist: %s', ksDir);
end

timeFile = fullfile(ksDir, 'spike_seconds_adj.npy');
unitFile = fullfile(ksDir, 'spike_clusters.npy');

if ~isfile(timeFile)
    error('Missing file: %s', timeFile);
end

if ~isfile(unitFile)
    error('Missing file: %s', unitFile);
end

spike_times_sec = readNPY(timeFile);
spike_clusters  = readNPY(unitFile);

spike_times_sec = double(spike_times_sec(:));
spike_clusters  = double(spike_clusters(:));

if numel(spike_times_sec) ~= numel(spike_clusters)
    error(['spike_seconds_adj.npy and spike_clusters.npy have different lengths ', ...
        'in %s'], ksDir);
end

spike_unit_time = [spike_clusters, spike_times_sec];
end


function eventunit_time = get_eventst_GLX(fn, sr_nidq)
% Extract rising-edge events from a NIDQ .bin file.
%
% Output:
%   eventunit_time(:,1) = 2000
%   eventunit_time(:,2) = event time in seconds
%
% Current assumptions:
% - the NIDQ file contains 4 channels total
% - the last channel is the packed digital channel
% - bit 2 is the event marker of interest
% - bit 1 would correspond to SpikeGLX/NIDQ sync, but is not used here

if exist(fn, 'file') ~= 2
    [filename, loadpath] = uigetfile('*.bin', 'Please select .bin file');
    if isequal(filename, 0)
        error('No nidq file was selected.');
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



counter = 0;
fseek(fid, 0, 'bof');

while ~feof(fid)
    d_size_t = min(d_size, floor((n_sample - counter) / n_chan));

    if d_size_t > 0
        [data, add_count] = fread(fid, [n_chan d_size_t], '*int16');

        % Use bit 2 from the packed digital channel as the event signal
        sync_channel(1, t:t+d_size_t-1) = bitget(data(end, :), 2);

        counter = counter + add_count;
        t = t + d_size_t;
    else
        break
    end
end

assert(counter == n_sample, 'Could not read all samples');

% Rising edges define event onsets
eventsample = find(diff(sync_channel) > 0.5);
eventt = eventsample / sr_nidq;
eventt = double(eventt(:));

n_events = numel(eventt);

fprintf('Assumed nidq sampling rate: %.10g\n', sr_nidq);
fprintf('Detected %d events\n', n_events);

eventunit_time = [ones(n_events, 1) * 2000, eventt];
end


function meta_file = find_nidq_meta_file(catgt_folder, destDir, run_g)
% Find a .nidq.meta file.
%
% Search order:
% 1) look directly inside catgt_folder
% 2) if none found, look directly inside destDir
%
% No recursive search is used.
% If multiple files are found, ask the user to choose one.

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


function sr_nidq = get_niSampRate(meta_file)
% Read niSampRate from a .meta file.

if ~isfile(meta_file)
    error('Meta file does not exist: %s', meta_file);
end

fid = fopen(meta_file, 'r');
if fid == -1
    error('Cannot open meta file: %s', meta_file);
end



sr_nidq = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));

    if ~ischar(line) || isempty(line)
        continue;
    end

    if startsWith(line, 'niSampRate=')
        val_str = strtrim(line(numel('niSampRate=')+1:end));
        sr_nidq = str2double(val_str);

        if isnan(sr_nidq)
            error('niSampRate in %s is not numeric: %s', meta_file, val_str);
        end
        return;
    end
end

error('niSampRate was not found in meta file: %s', meta_file);
end


function [expo_dur_time, nidq_dur_time] = run_exponidq_sanity_check(dur, expo_syncs, catgt_folder, sr_nidq)
% Compare Expo and NIDQ event timing within one block.
%
% Inputs:
%   dur        : pairwise event-to-event intervals from NIDQ, in seconds
%                usually computed as diff(event times) for one segmented block
%   expo_syncs : Expo sync/event times for the corresponding block
%   catgt_folder: output folder for the diagnostic plot
%
% Outputs:
%   expo_dur_time : Expo inter-event intervals in milliseconds
%   nidq_dur_time : NIDQ inter-event intervals in milliseconds
%
% What this function does:
% - converts Expo and NIDQ inter-event intervals to milliseconds
% - compares the up/down pattern of interval changes between Expo and NIDQ
% - warns if too many mismatches are detected
% - saves a diagnostic figure for visual inspection
%
% Notes:
% - Expo is assumed to use an internal timebase of 10000 Hz
% - NIDQ durations are passed in as seconds and converted here to ms
% - This is a sanity check for marker alignment, not a formal fitting procedure

sr_expo = 10000;

% Expo event intervals (in Expo ticks)
expo_dur = diff(double(expo_syncs));

% Convert both interval series to milliseconds
nidq_dur_time = dur .* 1000;
expo_dur_time = expo_dur ./ sr_expo .* 1000;

%% Check whether the change pattern matches between Expo and NIDQ
% If too many up/down transitions disagree, markers may be mismatched

threshold = 1;
numMismatchthre = 10;

dx = diff(expo_dur_time);
dy = diff(nidq_dur_time);

% Ignore tiny changes
dx(abs(dx) < threshold) = 0;
dy(abs(dy) < threshold) = 0;

sgnX = sign(dx);
sgnY = sign(dy)';

% Compare only locations where both sides show a nonzero change
valid = (sgnX ~= 0) & (sgnY ~= 0);

% Count sign mismatches
ind = find(valid & (sgnX ~= sgnY));
numMismatch = numel(ind);

fprintf(' number of mismatch for expo and nidq pulse events dur up down: %i\n', numMismatch);

if numMismatch > numMismatchthre
    warning('need check whether markers are correct');
end

%% Diagnostic plot
figure
plot(expo_dur_time)
hold on
plot(nidq_dur_time)
legend('expo', 'nidq')
saveas(gcf, fullfile(catgt_folder, ['events_check' num2str(numel(expo_syncs)) '.png']));

%% Display assumed timing information
fprintf('Assumed internal sampling rate expo: %i\n', sr_expo);
fprintf('Assumed sampling rate nidq: %i\n', sr_nidq);

end


function unit_run_metrics = plot_kilosort_run_metrics(spike_unit_time_trial, Allunits, stimTag, prestim_t, poststim_t, ksDir)
% Plot run-by-run summary metrics for one kilosort folder.
%
% Rows:
%   1) stimulus-period firing rate distribution
%   2) mean PSTH across unitsfr
%   3) d-prime distribution (stim vs prestim)
%   4) Fano factor distribution (classic FF, shortest stim duration)
%   5) Rsc distribution (classic spike-count correlation, shortest stim duration)
%
% Also returns unit_run_metrics:
%   unit_run_metrics{j}.unit_ids
%   unit_run_metrics{j}.fr_stim
%   unit_run_metrics{j}.dprime
%   unit_run_metrics{j}.fano_factor

Allunits = unique(Allunits(:));
Allunits = Allunits(Allunits ~= 2000);

nRun  = numel(spike_unit_time_trial);
nUnit = numel(Allunits);

unit_run_metrics = cell(nRun, 1);

binSize  = 0.001;   % 1 ms
gauss_sd = 0.010;   % 10 ms
gk = make_gaussian_kernel(binSize, gauss_sd);

% Fixed edges for easier visual comparison
fr_edges  = 0:2:60;
dp_edges  = -0.5:0.25:4;
ff_edges  = 0:0.25:5;
rsc_edges = -1:0.1:1;

figW = max(1200, 330 * nRun);
figH = 1400;
hfig = figure('Color', 'w', 'Position', [50 50 figW figH]);

for j = 1:nRun

    this_run_trials = spike_unit_time_trial{j};
    nTrial = numel(this_run_trials);

    % ------------------------------------------------------------
    % Extract stim end times from the second 2000 marker in each trial
    % ------------------------------------------------------------
    stim_end_all = nan(nTrial, 1);

    for m = 1:nTrial
        tr = this_run_trials{m};
        marker_t = tr(tr(:,1) == 2000, 2);

        
            stim_end_all(m) = marker_t(2);
      
    end

    valid_trial_idx = find(~isnan(stim_end_all) & stim_end_all > 0);

    if isempty(valid_trial_idx)

        unit_run_metrics{j} = struct();
        unit_run_metrics{j}.stim_tag = stimTag{j};
        unit_run_metrics{j}.unit_ids = Allunits(:);
        unit_run_metrics{j}.fr_stim = nan(nUnit, 1);
        unit_run_metrics{j}.dprime = nan(nUnit, 1);
        unit_run_metrics{j}.fano_factor = nan(nUnit, 1);
        unit_run_metrics{j}.valid_trial_idx = [];
        unit_run_metrics{j}.stim_end_median = NaN;
        unit_run_metrics{j}.stim_end_min = NaN;

        for rr = 1:5
            subplot(5, nRun, (rr-1)*nRun + j);
            axis off
            text(0.5, 0.5, 'No valid trials', 'HorizontalAlignment', 'center');
            if rr == 1
                title(stimTag{j}, 'Interpreter', 'none');
            end
        end
        continue;
    end

    stim_end_all = stim_end_all(valid_trial_idx);
    stim_end_med = median(stim_end_all);
    stim_end_min = min(stim_end_all);

    % ------------------------------------------------------------
    % PSTH window:
    % use the longest trial; shorter trials are padded with NaN
    % ------------------------------------------------------------
    trial_end_all = stim_end_all + poststim_t;
    psth_tmin = -prestim_t;
    psth_tmax = max(trial_end_all);

    psth_edges   = psth_tmin:binSize:(psth_tmax + binSize);
    psth_centers = psth_edges(1:end-1) + binSize/2;
    nBin = numel(psth_centers);

    % ------------------------------------------------------------
    % Allocate per-unit metrics
    % ------------------------------------------------------------
    fr_stim = nan(nUnit, 1);
    dprime  = nan(nUnit, 1);
    ff      = nan(nUnit, 1);

    unit_psth = nan(nUnit, nBin);

    % for classic Rsc: shortest stim window across trials
    stim_count_mat_rsc = zeros(nUnit, numel(valid_trial_idx));

    % Optional: save per-trial intermediate values too
    fr_stim_trial_mat          = nan(nUnit, numel(valid_trial_idx));
    dprime_stim_rate_trial_mat = nan(nUnit, numel(valid_trial_idx));
    pre_rate_trial_mat         = nan(nUnit, numel(valid_trial_idx));
    ff_count_trial_mat         = nan(nUnit, numel(valid_trial_idx));

    % ------------------------------------------------------------
    % Compute metrics unit by unit
    % ------------------------------------------------------------
    for u = 1:nUnit
        uid = Allunits(u);

        fr_rates_trial          = zeros(numel(valid_trial_idx), 1);   % full stim duration
        dprime_stim_rates_trial = nan(numel(valid_trial_idx), 1);     % first prestim_t of stim
        pre_rates_trial         = nan(numel(valid_trial_idx), 1);     % prestim window
        stim_counts_ff          = zeros(numel(valid_trial_idx), 1);

        psth_trials = nan(numel(valid_trial_idx), nBin);

        for tt = 1:numel(valid_trial_idx)
            tr = this_run_trials{valid_trial_idx(tt)};
            stim_end = stim_end_all(tt);
            trial_end = stim_end + poststim_t;

            % spike times for this unit in this trial
            spk_t = tr(tr(:,1) == uid, 2);

            % ----- firing rate for row 1: use the full stim duration of this trial -----
            stim_count_actual = sum(spk_t >= 0 & spk_t < stim_end);
            fr_rates_trial(tt) = stim_count_actual / stim_end;

            % ----- d-prime for row 3: use equal-duration windows -----
            dprime_win = prestim_t;

            if stim_end >= dprime_win
                stim_count_dp = sum(spk_t >= 0 & spk_t < dprime_win);
                dprime_stim_rates_trial(tt) = stim_count_dp / dprime_win;

                pre_count = sum(spk_t >= -prestim_t & spk_t < 0);
                pre_rates_trial(tt) = pre_count / prestim_t;
            else
                dprime_stim_rates_trial(tt) = NaN;
                pre_rates_trial(tt) = NaN;
            end

            % ----- classic FF and Rsc: use shortest stim window -----
            stim_counts_ff(tt) = sum(spk_t >= 0 & spk_t < stim_end_min);
            stim_count_mat_rsc(u, tt) = stim_counts_ff(tt);

            % ----- PSTH over longest window; pad beyond actual trial end with NaN -----
            r = nan(1, nBin);

            valid_bins = psth_centers <= trial_end;
            r(valid_bins) = 0;

            spk_keep = spk_t(spk_t >= psth_tmin & spk_t <= trial_end);
            c = histcounts(spk_keep, psth_edges);

            r(valid_bins) = c(valid_bins) / binSize;  % convert to sp/s
            r = smooth_with_nan(r, gk);

            psth_trials(tt, :) = r;
        end

        % save trial-level intermediate results
        fr_stim_trial_mat(u, :)          = fr_rates_trial(:)';
        dprime_stim_rate_trial_mat(u, :) = dprime_stim_rates_trial(:)';
        pre_rate_trial_mat(u, :)         = pre_rates_trial(:)';
        ff_count_trial_mat(u, :)         = stim_counts_ff(:)';

        % (a) FR: mean of per-trial rates
        fr_stim(u) = mean(fr_rates_trial);

        % (b) d-prime: invalid => NaN
        dprime(u) = calc_dprime_from_rates(dprime_stim_rates_trial, pre_rates_trial);

        % (c) classic FF on shortest stim window; invalid => NaN
        mu_ff = mean(stim_counts_ff);
        if mu_ff > 0
            ff(u) = var(stim_counts_ff, 0) / mu_ff;
        else
            ff(u) = NaN;
        end

        % (d) unit PSTH = average across trials, ignoring NaN padding
        unit_psth(u, :) = mean_ignore_nan(psth_trials, 1);
    end

    % ------------------------------------------------------------
    % Save per-unit metrics for this run block
    % ------------------------------------------------------------
    unit_run_metrics{j} = struct();
    unit_run_metrics{j}.stim_tag = stimTag{j};
    unit_run_metrics{j}.unit_ids = Allunits(:);
    unit_run_metrics{j}.fr_stim = fr_stim;
    unit_run_metrics{j}.dprime = dprime;
    unit_run_metrics{j}.fano_factor = ff;
    unit_run_metrics{j}.valid_trial_idx = valid_trial_idx(:);
    unit_run_metrics{j}.stim_end_median = stim_end_med;
    unit_run_metrics{j}.stim_end_min = stim_end_min;

    % Optional trial-level intermediate values
    unit_run_metrics{j}.fr_stim_trial = fr_stim_trial_mat;
    unit_run_metrics{j}.dprime_stim_rate_trial = dprime_stim_rate_trial_mat;
    unit_run_metrics{j}.prestim_rate_trial = pre_rate_trial_mat;
    unit_run_metrics{j}.ff_spikecount_trial = ff_count_trial_mat;

    % ------------------------------------------------------------
    % Classic Rsc from shortest-window spike counts
    % Exclude invalid pairs (zero variance or non-finite correlation)
    % ------------------------------------------------------------
    rsc_vals = [];

    if nUnit >= 2
        for u1 = 1:nUnit-1
            x = stim_count_mat_rsc(u1, :);

            for u2 = u1+1:nUnit
                y = stim_count_mat_rsc(u2, :);

                if std(x, 0, 2) == 0 || std(y, 0, 2) == 0
                    continue;
                end

                r = corr(x(:), y(:), 'rows', 'complete');

                if isfinite(r)
                    rsc_vals(end+1, 1) = r; %#ok<AGROW>
                end
            end
        end
    end

    % ------------------------------------------------------------
    % Population PSTH across units
    % ------------------------------------------------------------
    mean_psth = mean_ignore_nan(unit_psth, 1);
    std_psth  = std_ignore_nan(unit_psth, 0, 1);

    n_unit_valid_per_bin = sum(~isnan(unit_psth), 1);
    sem_psth = std_psth ./ sqrt(n_unit_valid_per_bin);
    sem_psth(n_unit_valid_per_bin <= 1) = NaN;

    % Keep only finite values for the distributions
    fr_plot  = fr_stim(isfinite(fr_stim));
    dp_plot  = dprime(isfinite(dprime));
    ff_plot  = ff(isfinite(ff));
    rsc_plot = rsc_vals(isfinite(rsc_vals));

    % ------------------------------------------------------------
    % Row 1: stimulus-period firing rate distribution
    % ------------------------------------------------------------
    subplot(5, nRun, j);
    histogram(fr_plot, fr_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.30 0.50 0.90], 'EdgeColor', [0.30 0.50 0.90]);
    title(stimTag{j}, 'Interpreter', 'none');
    xlabel('Firing rate during stim (sp/s)');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(fr_plot);
    xlim([fr_edges(1) fr_edges(end)]);

    % ------------------------------------------------------------
    % Row 2: mean PSTH across units
    % ------------------------------------------------------------
    subplot(5, nRun, nRun + j);
    hold on

    valid_p = isfinite(mean_psth) & isfinite(sem_psth);
    xfill = psth_centers(valid_p);
    ylo   = mean_psth(valid_p) - sem_psth(valid_p);
    yhi   = mean_psth(valid_p) + sem_psth(valid_p);

    fill([xfill fliplr(xfill)], [ylo fliplr(yhi)], ...
        [0.75 0.83 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    plot(psth_centers, mean_psth, 'b', 'LineWidth', 2);

    xline(0, 'k--', 'LineWidth', 1);
    xline(stim_end_med, 'k--', 'LineWidth', 1);

    xlabel('Time (s)');
    if j == 1
        ylabel('Response (sp/s)');
    end
    box off

    % ------------------------------------------------------------
    % Row 3: d-prime distribution
    % ------------------------------------------------------------
    subplot(5, nRun, 2*nRun + j);
    histogram(dp_plot, dp_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.25 0.60 0.90], 'EdgeColor', [0.25 0.60 0.90]);
    xlabel('d-prime');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(dp_plot);
    xlim([dp_edges(1) dp_edges(end)]);

    % ------------------------------------------------------------
    % Row 4: Fano factor distribution
    % ------------------------------------------------------------
    subplot(5, nRun, 3*nRun + j);
    histogram(ff_plot, ff_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.20 0.70 0.50], 'EdgeColor', [0.20 0.70 0.50]);
    xlabel('Fano factor');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(ff_plot);
    xlim([ff_edges(1) ff_edges(end)]);

    % ------------------------------------------------------------
    % Row 5: Rsc distribution
    % ------------------------------------------------------------
    subplot(5, nRun, 4*nRun + j);
    histogram(rsc_plot, rsc_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.15 0.65 0.70], 'EdgeColor', [0.15 0.65 0.70]);
    xlabel('R_{sc}');
    if j == 1
        ylabel('Proportion of pairs');
    end
    add_mu_median_text(rsc_plot);
    xlim([rsc_edges(1) rsc_edges(end)]);
end

sgtitle(strrep(ksDir, '\', '\\'), 'Interpreter', 'none');

saveas(hfig, fullfile(ksDir, 'run_metrics_summary.png'));
savefig(hfig, fullfile(ksDir, 'run_metrics_summary.fig'));
close(hfig);
end

function gk = make_gaussian_kernel(binSize, sigma)
% Create a normalized Gaussian kernel.
% sigma and binSize are both in seconds.

halfWidth = ceil(4 * sigma / binSize);
x = (-halfWidth:halfWidth) * binSize;
gk = exp(-(x.^2) / (2 * sigma^2));
gk = gk / sum(gk);
end


function y = smooth_with_nan(x, gk)
% Gaussian smoothing while respecting NaN padding.

valid = ~isnan(x);
x0 = x;
x0(~valid) = 0;

num = conv(x0, gk, 'same');
den = conv(double(valid), gk, 'same');

y = num ./ den;
y(den == 0) = NaN;
end


function dp = calc_dprime_from_rates(stim_rates, pre_rates)
% d-prime from trial-by-trial firing rates.
%
% dp = (mu_stim - mu_pre) / sqrt(0.5*(var_stim + var_pre))
%
% Invalid values return NaN.

stim_rates = stim_rates(isfinite(stim_rates));
pre_rates  = pre_rates(isfinite(pre_rates));

if isempty(stim_rates) || isempty(pre_rates)
    dp = NaN;
    return;
end

mu_stim = mean(stim_rates);
mu_pre  = mean(pre_rates);

var_stim = var(stim_rates, 0);
var_pre  = var(pre_rates, 0);

denom = sqrt(0.5 * (var_stim + var_pre));

if ~isfinite(denom) || denom <= 0
    dp = NaN;
else
    dp = (mu_stim - mu_pre) / denom;
end
end


function m = mean_ignore_nan(X, dim)
% Mean ignoring NaNs.
% dim = 1 or 2

if nargin < 2
    dim = 1;
end

valid = ~isnan(X);
X0 = X;
X0(~valid) = 0;

n = sum(valid, dim);
m = sum(X0, dim) ./ n;
m(n == 0) = NaN;
end


function s = std_ignore_nan(X, flag, dim)
% Std ignoring NaNs.
% Supports dim = 1 or 2.

if nargin < 2 || isempty(flag)
    flag = 0;
end
if nargin < 3
    dim = 1;
end

m = mean_ignore_nan(X, dim);

if dim == 1
    mrep = repmat(m, size(X,1), 1);
else
    mrep = repmat(m, 1, size(X,2));
end

valid = ~isnan(X);
D = X - mrep;
D(~valid) = 0;

n = sum(valid, dim);

if flag == 0
    denom = n - 1;
else
    denom = n;
end

s = sqrt(sum(D.^2, dim) ./ denom);
s(denom <= 0) = NaN;
end


function add_mu_median_text(x)
% Add mean / median text to current axis.

x = x(isfinite(x));

if isempty(x)
    text(0.65, 0.90, 'n=0', 'Units', 'normalized', ...
        'FontSize', 10, 'VerticalAlignment', 'top');
    return;
end

mu = mean(x);
md = median(x);

text(0.62, 0.92, sprintf('n=%d\n\\mu=%.3f\nm=%.3f', numel(x), mu, md), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
end