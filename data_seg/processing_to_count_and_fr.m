%% =========================================================================
% processing_to_count_and_fr
%
% Purpose:
%   For selected run(s) identified by stim_tag, generate binned spike-count
%   and firing-rate data from spike_unit_time_trial, then compute several
%   normalized versions of the same data.
%
% Inputs required in each kilosort folder:
%   - spike_unit_time_trial.mat
%   - unit_condition_metrics.mat
%
% Main output saved in each kilosort folder:
%   - bined_data_allruns.mat
%
% Output structure:
%   bined_data_allruns{r}
%       .stim_tag
%       .unit_ids
%       .analysis_window
%       .bin_size
%       .bin_edges
%       .bin_centers
%       .condition_fields
%       .condition_index_per_trial
%       .conditions
%       .raw_count
%       .raw_fr
%       .z_within_trial
%       .z_within_condition
%       .z_across_conditions
%       .demean_count_within_trial
%       .demean_fr_within_trial
%       .demean_pooledsd_within_condition
%
% Notes:
%   1. All major data arrays are stored as:
%          unit x trial x bin
%   2. Different selected runs may use different analysis windows and
%      different bin sizes.
%   3. When normalization is undefined because the relevant standard
%      deviation is zero, this script returns NaN so that these cases can
%      be excluded later or handled explicitly in downstream analysis.
% =========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

root_folder = 'I:\np_data';
runName     = 'RafiL001p0120';
runind      = 1;          % run index after -g
probes      = [0,1];      % probe indices after -prb

% -------------------------------------------------------------------------
% One entry per requested run.
% Each requested stim_tag can have its own analysis_window and bin_size.
% -------------------------------------------------------------------------
run_specs = struct([]);

run_specs(1).stim_tag = '[RFG_coarse2dg_99_4_150isi]';
run_specs(1).analysis_window = [0 0.15];
run_specs(1).bin_size = 0.15;

run_specs(2).stim_tag = '[dir12_gpl_2_200isi_fixedphase]';
run_specs(2).analysis_window = [0 0.2];
run_specs(2).bin_size = 0.2;

run_specs(3).stim_tag = '_2[Gpl2_2c_2sz_400_2_200isi]';
run_specs(3).analysis_window = [0 0.4];
run_specs(3).bin_size = 0.02;

%% ----------------------- Build shared session paths -----------------------

run_g   = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);

fprintf('destDir : %s\n', destDir);

%% ----------------------- Process each probe folder -----------------------

for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr   = sprintf('imec%d', thisProbe);

    probe_folder = fullfile(destDir, ['catgt_' run_g], [run_g '_' imecStr]);

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    if ~isfolder(probe_folder)
        warning('probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    %% ----------------------- Find kilosort folders -----------------------

    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);

    if isempty(d)
        warning('No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf('Found %d kilosort folder(s) under probe %d.\n', numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------

    for i = 1:numel(d)

        ksDir = fullfile(d(i).folder, d(i).name);

        fprintf('\nProcessing probe %d, ksDir: %s\n', thisProbe, ksDir);

        try
            %% ----------------------- Load required files -----------------------

            trial_file = fullfile(ksDir, 'spike_unit_time_trial.mat');
            if ~isfile(trial_file)
                error('Missing file: %s', trial_file);
            end

            S = load(trial_file, 'spike_unit_time_trial', 'prestim_t', 'poststim_t');

            if ~isfield(S, 'spike_unit_time_trial')
                error('spike_unit_time_trial not found in %s', trial_file);
            end
            if ~isfield(S, 'prestim_t')
                error('prestim_t not found in %s', trial_file);
            end
            if ~isfield(S, 'poststim_t')
                error('poststim_t not found in %s', trial_file);
            end

            spike_unit_time_trial = S.spike_unit_time_trial;
            prestim_t = S.prestim_t;
            poststim_t = S.poststim_t;

            cond_file = fullfile(ksDir, 'unit_condition_metrics.mat');
            if ~isfile(cond_file)
                error('Missing file: %s', cond_file);
            end

            Sc = load(cond_file, 'unit_condition_metrics');
            if ~isfield(Sc, 'unit_condition_metrics')
                error('unit_condition_metrics not found in %s', cond_file);
            end

            unit_condition_metrics = Sc.unit_condition_metrics;

            if numel(spike_unit_time_trial) ~= numel(unit_condition_metrics)
                error(['Number of runs mismatch between spike_unit_time_trial (%d) ' ...
                    'and unit_condition_metrics (%d) in ksDir %s'], ...
                    numel(spike_unit_time_trial), numel(unit_condition_metrics), ksDir);
            end

            %% ----------------------- Process requested runs -----------------------

            bined_data_allruns = cell(numel(run_specs), 1);
            all_run_tags = get_all_run_tags(unit_condition_metrics);

            for r = 1:numel(run_specs)

                wanted_tag = run_specs(r).stim_tag;
                analysis_window = run_specs(r).analysis_window;
                bin_size = run_specs(r).bin_size;

                run_idx = find(strcmp(all_run_tags, wanted_tag));

                if isempty(run_idx)
                    error('Requested stim_tag not found in ksDir %s: %s', ksDir, wanted_tag);
                end

                if numel(run_idx) > 1
                    error('Duplicate stim_tag found in unit_condition_metrics for ksDir %s: %s', ...
                        ksDir, wanted_tag);
                end

                this_run_trials = spike_unit_time_trial{run_idx};
                this_run_metrics = unit_condition_metrics{run_idx};

                bined_data_allruns{r} = process_one_run( ...
                    this_run_trials, this_run_metrics, ...
                    analysis_window, bin_size, ...
                    prestim_t, poststim_t, run_idx, wanted_tag);
            end

            %% ----------------------- Save output -----------------------

            save(fullfile(ksDir, 'bined_data_allruns.mat'), ...
                'bined_data_allruns');

            fprintf('Saved:\n');
            fprintf('  %s\n', fullfile(ksDir, 'bined_data_allruns.mat'));

        catch ME
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', ME.message);
        end
    end
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================

function out = process_one_run(this_run_trials, this_run_metrics, analysis_window, bin_size, prestim_t, poststim_t, run_idx, stim_tag)
%% =========================================================================
% process_one_run
%
% Purpose:
%   For one selected run, convert trial-level spike times into binned count
%   and firing-rate matrices, then compute several normalized versions of
%   those matrices.
% =========================================================================

if ~isfield(this_run_metrics, 'unit_ids')
    error('unit_ids missing in unit_condition_metrics{%d}.', run_idx);
end
if ~isfield(this_run_metrics, 'condition_index_per_trial')
    error('condition_index_per_trial missing in unit_condition_metrics{%d}.', run_idx);
end
if ~isfield(this_run_metrics, 'conditions')
    error('conditions missing in unit_condition_metrics{%d}.', run_idx);
end
if ~isfield(this_run_metrics, 'condition_fields')
    error('condition_fields missing in unit_condition_metrics{%d}.', run_idx);
end

unit_ids = this_run_metrics.unit_ids(:);
condition_index_per_trial = this_run_metrics.condition_index_per_trial(:);
condition_fields = this_run_metrics.condition_fields(:);

conditions_identity = extract_condition_identity(this_run_metrics.conditions, condition_fields);

nUnit = numel(unit_ids);
nTrial = numel(this_run_trials);

if numel(condition_index_per_trial) ~= nTrial
    error(['condition_index_per_trial length (%d) does not match number of trials (%d) ' ...
        'for run %s'], numel(condition_index_per_trial), nTrial, stim_tag);
end

[nBin, bin_edges, bin_centers] = make_bin_definition(analysis_window, bin_size);

stim_end_all = nan(nTrial, 1);
for t = 1:nTrial
    stim_end_all(t) = get_trial_stim_end(this_run_trials{t});
end

validate_analysis_window(analysis_window, prestim_t, poststim_t, stim_end_all, stim_tag);

% -------------------------------------------------------------------------
% Build raw count matrix: unit x trial x bin
% -------------------------------------------------------------------------
raw_count = zeros(nUnit, nTrial, nBin);

for t = 1:nTrial
    tr = this_run_trials{t};

    for u = 1:nUnit
        uid = unit_ids(u);
        spk_t = tr(tr(:,1) == uid, 2);

        if isempty(spk_t)
            c = zeros(1, nBin);
        else
            spk_keep = spk_t(spk_t >= analysis_window(1) & spk_t <= analysis_window(2));
            c = histcounts(spk_keep, bin_edges);
        end

        raw_count(u, t, :) = reshape(c, [1 1 nBin]);
    end
end

raw_count = double(raw_count);
raw_fr = raw_count ./ bin_size;

% -------------------------------------------------------------------------
% Normalized data
% -------------------------------------------------------------------------

z_within_trial = zscore_within_trial(raw_count);
z_within_condition = zscore_within_condition(raw_count, condition_index_per_trial);
z_across_conditions = zscore_across_conditions(raw_count);

demean_count_within_trial = demean_within_trial(raw_count);
demean_fr_within_trial    = demean_within_trial(raw_fr);

demean_pooledsd_within_condition_data = demean_pooledsd_within_condition(raw_count, condition_index_per_trial);

% -------------------------------------------------------------------------
% Pack output
% -------------------------------------------------------------------------
out = struct();
out.stim_tag = stim_tag;

out.unit_ids = unit_ids;
out.analysis_window = analysis_window;
out.bin_size = bin_size;
out.bin_edges = bin_edges;
out.bin_centers = bin_centers;

out.condition_fields = condition_fields;
out.condition_index_per_trial = condition_index_per_trial;
out.conditions = conditions_identity;

out.raw_count = raw_count;
out.raw_fr = raw_fr;

out.z_within_trial = z_within_trial;
out.z_within_condition = z_within_condition;
out.z_across_conditions = z_across_conditions;

out.demean_count_within_trial = demean_count_within_trial;
out.demean_fr_within_trial = demean_fr_within_trial;

out.demean_pooledsd_within_condition = demean_pooledsd_within_condition_data;

end


function [nBin, bin_edges, bin_centers] = make_bin_definition(analysis_window, bin_size)
t0 = analysis_window(1);
t1 = analysis_window(2);
span = t1 - t0;

nBin_float = span / bin_size;
nBin = round(nBin_float);

tol = 1e-10;
if abs(nBin_float - nBin) > tol
    error(['analysis_window span (%.12g s) is not an integer multiple of bin_size ' ...
        '(%.12g s).'], span, bin_size);
end

bin_edges = t0 + (0:nBin) * bin_size;
bin_edges(end) = t1;
bin_centers = bin_edges(1:end-1) + bin_size/2;
end


function validate_analysis_window(analysis_window, prestim_t, poststim_t, stim_end_all, stim_tag)

if analysis_window(1) < -prestim_t
    error(['analysis_window start %.6g is earlier than available prestim coverage ' ...
        '(-prestim_t = %.6g) for run %s.'], ...
        analysis_window(1), -prestim_t, stim_tag);
end

latest_safe_end = min(stim_end_all + poststim_t);

if analysis_window(2) > latest_safe_end
    error(['analysis_window end %.6g exceeds the shortest available trial coverage ' ...
        '(min(stim_end + poststim_t) = %.6g) for run %s.'], ...
        analysis_window(2), latest_safe_end, stim_tag);
end

end


function stim_end = get_trial_stim_end(tr)

if isempty(tr)
    error('Encountered an empty trial matrix.');
end

marker_t = tr(tr(:,1) == 2000, 2);

if numel(marker_t) < 2
    error('A trial does not contain at least two event markers with unit ID 2000.');
end

stim_end = marker_t(2);
end


function all_tags = get_all_run_tags(unit_condition_metrics)

all_tags = cell(numel(unit_condition_metrics), 1);

for j = 1:numel(unit_condition_metrics)
    if ~isfield(unit_condition_metrics{j}, 'stim_tag')
        error('stim_tag missing in unit_condition_metrics{%d}.', j);
    end
    all_tags{j} = unit_condition_metrics{j}.stim_tag;
end

end


function conditions_identity = extract_condition_identity(conditions_in, condition_fields)

nCond = numel(conditions_in);
conditions_identity = repmat(struct(), nCond, 1);

for c = 1:nCond
    conditions_identity(c).trial_indices = conditions_in(c).trial_indices(:);

    for k = 1:numel(condition_fields)
        f = condition_fields{k};
        if isfield(conditions_in(c), f)
            conditions_identity(c).(f) = conditions_in(c).(f);
        else
            error('Condition field %s is missing in conditions(%d).', f, c);
        end
    end
end

end


function Z = zscore_within_trial(X)

[nUnit, nTrial, nBin] = size(X);
Z = zeros(nUnit, nTrial, nBin);

for u = 1:nUnit
    for t = 1:nTrial
        v = reshape(X(u, t, :), [nBin, 1]);
        Z(u, t, :) = reshape(zscore_vector_safe(v), [1 1 nBin]);
    end
end

end


function Z = zscore_within_condition(X, condition_index_per_trial)

[nUnit, ~, nBin] = size(X);
Z = zeros(size(X));

cond_ids = unique(condition_index_per_trial(:))';
cond_ids = cond_ids(isfinite(cond_ids));

for u = 1:nUnit
    for c = cond_ids
        trial_idx = find(condition_index_per_trial == c);

        pooled = reshape(X(u, trial_idx, :), [], 1);
        z_pooled = zscore_vector_safe(pooled);
        Z(u, trial_idx, :) = reshape(z_pooled, [1 numel(trial_idx) nBin]);
    end
end

end


function Z = zscore_across_conditions(X)

[nUnit, nTrial, nBin] = size(X);
Z = zeros(nUnit, nTrial, nBin);

for u = 1:nUnit
    pooled = reshape(X(u, :, :), [], 1);
    z_pooled = zscore_vector_safe(pooled);
    Z(u, :, :) = reshape(z_pooled, [1 nTrial nBin]);
end

end


function Y = demean_within_trial(X)

[nUnit, nTrial, nBin] = size(X);
Y = zeros(nUnit, nTrial, nBin);

for u = 1:nUnit
    for t = 1:nTrial
        v = reshape(X(u, t, :), [nBin, 1]);
        v = v - mean(v);
        Y(u, t, :) = reshape(v, [1 1 nBin]);
    end
end

end


function Y = demean_pooledsd_within_condition(X, condition_index_per_trial)
%% =========================================================================
% demean_pooledsd_within_condition
%
% Purpose:
%   For each neuron and each condition:
%     1) subtract the mean across bins within each trial
%     2) pool all demeaned points across trials and bins in that condition
%     3) divide by the pooled standard deviation
%
% Input:
%   X : numeric array
%       Data array of size unit x trial x bin.
%   condition_index_per_trial : numeric vector
%       Condition label for each trial.
%
% Output:
%   Y : numeric array
%       Output array of the same size as X.
%
% Notes:
%   If the pooled standard deviation is zero, the output for that neuron
%   and condition is returned as NaN.
% =========================================================================

[nUnit, ~, nBin] = size(X);
Y = zeros(size(X));

cond_ids = unique(condition_index_per_trial(:))';
cond_ids = cond_ids(isfinite(cond_ids));

for u = 1:nUnit
    for c = cond_ids
        trial_idx = find(condition_index_per_trial == c);
        nThisTrial = numel(trial_idx);

        tmp = nan(nThisTrial, nBin);

        for tt = 1:nThisTrial
            tr = trial_idx(tt);
            v = reshape(X(u, tr, :), [nBin, 1]);
            v = v - mean(v);
            tmp(tt, :) = v(:)';
        end

        pooled = tmp(:);
        sd_pooled = std(pooled, 0);

        if sd_pooled == 0
            tmp_out = nan(size(tmp));
        else
            tmp_out = tmp ./ sd_pooled;
        end

        Y(u, trial_idx, :) = reshape(tmp_out, [1 nThisTrial nBin]);
    end
end

end


function z = zscore_vector_safe(v)

mu = mean(v);
sd = std(v, 0);

if sd == 0
    z = NaN(size(v));
else
    z = (v - mu) ./ sd;
end

end