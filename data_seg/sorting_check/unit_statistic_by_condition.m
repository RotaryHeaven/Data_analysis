%% =========================================================================
% Condition-based run metrics from spike_unit_time_trial + stiminfo
%
% Purpose:
%   For each run (stimulus block), group trials into stimulus conditions
%   using stiminfo (ignoring condition_ids), then compute condition-level
%   metrics for each unit:
%       - firing rate during stimulus
%       - d-prime (stim vs prestim)
%       - Fano factor
%
%   Finally, for plotting, take across-condition summary values for each
%   unit within each run:
%       - maximum FR
%       - minimum FR
%       - maximum d-prime
%       - maximum Fano factor
%
% Inputs required:
%   1. output in each kilosort folder:
%        - spike_unit_time_trial.mat
%        - unit_run_metrics.mat
%   2. output in catgt_folder:
%        - stiminfo.mat
%
% Expected stiminfo format:
%   stiminfo{1:nRun}   = per-run stimulus info structs
%   stiminfo{nRun+1}   = stimTag (cell/string array of run names)
%
% Outputs saved in each kilosort folder:
%   - unit_condition_metrics.mat
%   - condition_metrics_summary_max.png
%   - condition_metrics_summary_max.fig
%
% =========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

root_folder = 'I:\np_data';
runName     = 'RafiL001p0120';
runind      = 1;          % run index after -g
probes      = [0,1];      % probe indices after -prb

%% ----------------------- Build shared session paths -----------------------

% Build the session name with the g-index suffix.
run_g   = sprintf('%s_g%d', runName, runind);
destDir = fullfile(root_folder, run_g);

% Define CatGT folder and stiminfo file path.
catgt_folder  = fullfile(destDir, ['catgt_' run_g]);
stiminfo_file = fullfile(catgt_folder, 'stiminfo.mat');

fprintf('destDir      : %s\n', destDir);
fprintf('catgt_folder : %s\n', catgt_folder);
fprintf('stiminfo_file: %s\n', stiminfo_file);

% stiminfo.mat must exist before processing.
if ~isfile(stiminfo_file)
    error('stiminfo.mat does not exist: %s', stiminfo_file);
end

% Load only the stiminfo variable.
Sstim = load(stiminfo_file, 'stiminfo');
if ~isfield(Sstim, 'stiminfo')
    error('stiminfo not found in %s', stiminfo_file);
end

stiminfo_all = Sstim.stiminfo;

% Expected format:
%   stiminfo{1:end-1} = run-wise stiminfo structs
%   stiminfo{end}     = stimTag
if ~iscell(stiminfo_all) || numel(stiminfo_all) < 2
    error(['stiminfo.mat must contain a cell array named stiminfo with ', ...
        'format stiminfo{1:nRun}=structs and stiminfo{nRun+1}=stimTag.']);
end

stiminfo = stiminfo_all(1:end-1);
stimTag_from_stiminfo = stiminfo_all{end};

% stimTag must be stored as cell array or string array.
if ~iscell(stimTag_from_stiminfo) && ~isstring(stimTag_from_stiminfo)
    error(['stiminfo{end} must be stimTag, stored as a cell array or ', ...
        'string array of run names.']);
end

% Normalize stimTag into a standard cell array of char vectors.
stimTag_from_stiminfo = normalize_stimtag_list(stimTag_from_stiminfo);

% The number of run structs must match the number of run names.
if numel(stimTag_from_stiminfo) ~= numel(stiminfo)
    error(['stiminfo format error: number of run structs is %d, but ', ...
        'number of run names in stiminfo{end} is %d.'], ...
        numel(stiminfo), numel(stimTag_from_stiminfo));
end

% Each run entry must be a struct.
for j = 1:numel(stiminfo)
    if ~isstruct(stiminfo{j})
        error('stiminfo{%d} is not a struct.', j);
    end
end

%% ----------------------- Process each probe folder -----------------------

for ip = 1:numel(probes)

    thisProbe = probes(ip);
    imecStr   = sprintf('imec%d', thisProbe);

    % Build the current probe folder path.
    probe_folder = fullfile( ...
        destDir, ...
        ['catgt_' run_g], ...
        [run_g '_' imecStr] ...
        );

    fprintf('\n============================================================\n');
    fprintf('Processing probe %d\n', thisProbe);
    fprintf('probe_folder: %s\n', probe_folder);
    fprintf('============================================================\n');

    % Skip this probe if the folder does not exist.
    if ~isfolder(probe_folder)
        warning('probe_folder does not exist, skipping probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    %% ----------------------- Find all kilosort folders -----------------------

    % Search for all kilosort* directories under the probe folder.
    d = dir(fullfile(probe_folder, 'kilosort*'));
    d = d([d.isdir]);

    if isempty(d)
        warning('No kilosort* folders found under probe %d: %s', ...
            thisProbe, probe_folder);
        continue;
    end

    % Sort folder names alphabetically for reproducible processing order.
    [~, idx] = sort(lower({d.name}));
    d = d(idx);

    fprintf('Found %d kilosort folder(s) under probe %d.\n', numel(d), thisProbe);

    %% ----------------------- Process each kilosort folder -----------------------

    for i = 1:numel(d)

        ksDir = fullfile(d(i).folder, d(i).name);

        % If all outputs already exist, skip this folder.
        done_file_1 = fullfile(ksDir, 'unit_condition_metrics.mat');
        done_file_2 = fullfile(ksDir, 'condition_metrics_summary_maxmin.png');
        done_file_3 = fullfile(ksDir, 'condition_metrics_summary_maxmin.fig');

        if isfile(done_file_1) && isfile(done_file_2) && isfile(done_file_3)
            fprintf('\nSkipping ksDir (already analyzed): %s\n', ksDir);
            continue;
        end

        fprintf('\nProcessing probe %d, ksDir: %s\n', thisProbe, ksDir);

        try
            %% ----------------------- Load spike_unit_time_trial -----------------------

            % Load trial-aligned spike data and stimulus timing parameters.
            trial_file = fullfile(ksDir, 'spike_unit_time_trial.mat');
            if ~isfile(trial_file)
                error('Missing file: %s', trial_file);
            end

            S = load(trial_file, 'spike_unit_time_trial', 'prestim_t', 'poststim_t');

            if ~isfield(S, 'spike_unit_time_trial')
                error('spike_unit_time_trial not found in %s', trial_file);
            end

            spike_unit_time_trial = S.spike_unit_time_trial;

            if ~isfield(S, 'prestim_t')
                error('prestim_t not found in %s', trial_file);
            end
            prestim_t = S.prestim_t;

            % poststim_t is optional.
            if isfield(S, 'poststim_t')
                poststim_t = S.poststim_t;
            else
                poststim_t = NaN;
            end

            %% ----------------------- Load unit_run_metrics from Program 1 -----------------------

            % Load Program 1 output for run-wise unit information.
            run_metrics_file = fullfile(ksDir, 'unit_run_metrics.mat');
            if ~isfile(run_metrics_file)
                error('Missing file: %s', run_metrics_file);
            end

            Srun = load(run_metrics_file, 'unit_run_metrics');
            if ~isfield(Srun, 'unit_run_metrics')
                error('unit_run_metrics not found in %s', run_metrics_file);
            end

            unit_run_metrics = Srun.unit_run_metrics;

            %% ----------------------- Cross-check run counts -----------------------

            % The number of runs must match across spike data, run metrics,
            % and stiminfo.
            if numel(spike_unit_time_trial) ~= numel(stiminfo)
                error(['Number of runs mismatch between spike_unit_time_trial (%d) ' ...
                    'and stiminfo (%d) in ksDir %s'], ...
                    numel(spike_unit_time_trial), numel(stiminfo), ksDir);
            end

            if numel(unit_run_metrics) ~= numel(stiminfo)
                error(['Number of runs mismatch between unit_run_metrics (%d) ' ...
                    'and stiminfo (%d) in ksDir %s'], ...
                    numel(unit_run_metrics), numel(stiminfo), ksDir);
            end

            %% ----------------------- Cross-check run names one by one -----------------------

            % Verify that Program 1 and Program 2 use the same run labels.
            for j = 1:numel(stiminfo)

                if ~isfield(unit_run_metrics{j}, 'stim_tag')
                    error('unit_run_metrics{%d}.stim_tag is missing in %s', j, run_metrics_file);
                end

                stim_tag_from_prog1 = normalize_one_stimtag(unit_run_metrics{j}.stim_tag);
                stim_tag_from_prog2 = stimTag_from_stiminfo{j};

                if ~strcmp(stim_tag_from_prog1, stim_tag_from_prog2)
                    error(['Run-name mismatch at run %d in ksDir %s\n' ...
                        'Program 1 unit_run_metrics{%d}.stim_tag = %s\n' ...
                        'Program 2 stiminfo{end}{%d}            = %s'], ...
                        j, ksDir, j, stim_tag_from_prog1, j, stim_tag_from_prog2);
                end
            end

            fprintf('Run count and run-name checks passed for ksDir.\n');

            %% ----------------------- Compute condition-level metrics -----------------------

            % One output struct per run.
            unit_condition_metrics = cell(numel(spike_unit_time_trial), 1);

            for j = 1:numel(spike_unit_time_trial)

                if ~isfield(unit_run_metrics{j}, 'unit_ids')
                    error('unit_run_metrics{%d}.unit_ids is missing in %s', j, run_metrics_file);
                end

                % Get unit IDs and run label for the current run.
                Allunits_thisrun = unit_run_metrics{j}.unit_ids(:);
                stim_tag_thisrun = stimTag_from_stiminfo{j};

                % Compute condition-level metrics for this run.
                unit_condition_metrics{j} = compute_run_condition_metrics( ...
                    spike_unit_time_trial{j}, stiminfo{j}, ...
                    Allunits_thisrun, prestim_t, poststim_t, ...
                    j, stim_tag_thisrun);
            end

            % Save computed condition-level metrics.
            save(fullfile(ksDir, 'unit_condition_metrics.mat'), ...
                'unit_condition_metrics', 'prestim_t', 'poststim_t');

            fprintf('Saved:\n');
            fprintf('  %s\n', fullfile(ksDir, 'unit_condition_metrics.mat'));

            %% ----------------------- Plot summary using max/min across conditions -----------------------

            % Plot per-run summary histograms using across-condition
            % summary values for each unit.
            plot_condition_summary_maxmin(unit_condition_metrics, ksDir);

            fprintf('  %s\n', fullfile(ksDir, 'condition_metrics_summary_maxmin.png'));
            fprintf('  %s\n', fullfile(ksDir, 'condition_metrics_summary_maxmin.fig'));

        catch ME
            % Report error but continue with the next ks folder.
            fprintf(2, 'Error in probe %d, ksDir %s\n', thisProbe, ksDir);
            fprintf(2, '%s\n', ME.message);
        end
    end
end

fprintf('\nDone.\n');

%% ======================= Local functions =======================


function run_metrics = compute_run_condition_metrics(this_run_trials, stiminfo_run, Allunits, prestim_t, poststim_t, run_idx, stim_tag)
%% =========================================================================
% compute_run_condition_metrics
%
% Purpose:
%   For one run, group trials into stimulus conditions and compute
%   condition-level metrics for each unit, including:
%       - firing rate during the full stimulus window
%       - d-prime between stimulus and prestimulus activity
%       - Fano factor using a common shortest stimulus window
%
% Inputs:
%   this_run_trials : cell array
%       Trial-wise spike/event matrix for one run.
%   stiminfo_run : struct
%       Trial-wise stimulus information for one run.
%   Allunits : numeric vector
%       Unit IDs to be analyzed in this run.
%   prestim_t : numeric scalar
%       Prestimulus duration used for d-prime calculation.
%   poststim_t : numeric scalar
%       Poststimulus duration loaded from file (stored in output).
%   run_idx : numeric scalar
%       Run index.
%   stim_tag : char/string
%       Run label.
%
% Outputs:
%   run_metrics : struct
%       Structure containing per-condition metrics, trial-level intermediates,
%       and per-unit summary values across conditions.
%
% =========================================================================

% Number of trials from spike data and from stiminfo.
nTrial_spk  = numel(this_run_trials);
nTrial_stim = get_n_trial_from_stiminfo(stiminfo_run);

% Trial counts must match.
if nTrial_spk ~= nTrial_stim
    error(['Trial count mismatch in run %d: spike_unit_time_trial has %d trials, ' ...
        'stiminfo has %d trials.'], run_idx, nTrial_spk, nTrial_stim);
end

% Extract stimulus end time for each trial.
stim_end_all = nan(nTrial_spk, 1);
for m = 1:nTrial_spk
    stim_end_all(m) = get_trial_stim_end(this_run_trials{m});
end

% Group trials into conditions using trial-wise stiminfo fields while
% ignoring condition_ids.
[condIndex, condDefs, used_fields] = group_trials_by_condition(stiminfo_run);

nCond = numel(condDefs);
nUnit = numel(Allunits);

% Preallocate unit-by-condition summary matrices.
fr_by_cond = nan(nUnit, nCond);
dp_by_cond = nan(nUnit, nCond);
ff_by_cond = nan(nUnit, nCond);

% Initialize output structure.
run_metrics = struct();
run_metrics.run_index = run_idx;
run_metrics.stim_tag = stim_tag;
run_metrics.n_trials = nTrial_spk;
run_metrics.unit_ids = Allunits(:);
run_metrics.condition_fields = used_fields(:);
run_metrics.condition_index_per_trial = condIndex(:);
run_metrics.n_conditions = nCond;

run_metrics.conditions = condDefs;

run_metrics.prestim_t = prestim_t;
run_metrics.poststim_t = poststim_t;

for c = 1:nCond

    % Trial indices belonging to this condition.
    trial_idx = condDefs(c).trial_indices(:);

    % Defensive handling for empty conditions.
    if isempty(trial_idx)
        run_metrics.conditions(c).stim_end_min = NaN;
        run_metrics.conditions(c).stim_end_median = NaN;
        run_metrics.conditions(c).fr_stim = nan(nUnit, 1);
        run_metrics.conditions(c).dprime = nan(nUnit, 1);
        run_metrics.conditions(c).fano_factor = nan(nUnit, 1);
        run_metrics.conditions(c).fr_stim_trial = nan(nUnit, 0);
        run_metrics.conditions(c).dprime_stim_rate_trial = nan(nUnit, 0);
        run_metrics.conditions(c).prestim_rate_trial = nan(nUnit, 0);
        run_metrics.conditions(c).ff_spikecount_trial = nan(nUnit, 0);
        continue;
    end

    % Use condition-specific stimulus duration summaries.
    stim_end_cond = stim_end_all(trial_idx);
    stim_end_min  = min(stim_end_cond);
    stim_end_med  = median(stim_end_cond);

    % Preallocate per-unit outputs for this condition.
    fr_stim = nan(nUnit, 1);
    dprime  = nan(nUnit, 1);
    ff      = nan(nUnit, 1);

    % Store trial-level intermediate results for later inspection.
    fr_stim_trial_mat          = nan(nUnit, numel(trial_idx));
    dprime_stim_rate_trial_mat = nan(nUnit, numel(trial_idx));
    pre_rate_trial_mat         = nan(nUnit, numel(trial_idx));
    ff_count_trial_mat         = nan(nUnit, numel(trial_idx));

    for u = 1:nUnit
        uid = Allunits(u);

        % Temporary per-trial vectors for one unit.
        fr_rates_trial          = nan(numel(trial_idx), 1);
        dprime_stim_rates_trial = nan(numel(trial_idx), 1);
        pre_rates_trial         = nan(numel(trial_idx), 1);
        stim_counts_ff          = nan(numel(trial_idx), 1);

        for tt = 1:numel(trial_idx)
            tr = this_run_trials{trial_idx(tt)};
            stim_end = stim_end_all(trial_idx(tt));

            % Extract spike times for the current unit in this trial.
            spk_t = tr(tr(:,1) == uid, 2);

            % (1) Firing rate over the full stimulus window of this trial.
            stim_count_actual = sum(spk_t >= 0 & spk_t < stim_end);
            fr_rates_trial(tt) = stim_count_actual / stim_end;

            % (2) d-prime uses matched prestim/stim windows of length prestim_t.
            if stim_end >= prestim_t
                stim_count_dp = sum(spk_t >= 0 & spk_t < prestim_t);
                dprime_stim_rates_trial(tt) = stim_count_dp / prestim_t;

                pre_count = sum(spk_t >= -prestim_t & spk_t < 0);
                pre_rates_trial(tt) = pre_count / prestim_t;
            end

            % (3) Fano factor uses the shortest stimulus window within the
            % current condition so that spike counts are comparable.
            stim_counts_ff(tt) = sum(spk_t >= 0 & spk_t < stim_end_min);
        end

        % Save trial-level matrices.
        fr_stim_trial_mat(u, :)          = fr_rates_trial(:)';
        dprime_stim_rate_trial_mat(u, :) = dprime_stim_rates_trial(:)';
        pre_rate_trial_mat(u, :)         = pre_rates_trial(:)';
        ff_count_trial_mat(u, :)         = stim_counts_ff(:)';

        % Compute final unit-level condition metrics.
        % Compute final unit-level condition metrics WITHOUT ignoring NaN.
        fr_stim(u) = mean(fr_rates_trial);
        dprime(u)  = calc_dprime_from_rates(dprime_stim_rates_trial, pre_rates_trial);

        mu_ff = mean(stim_counts_ff);
        if isfinite(mu_ff) && mu_ff > 0
            ff(u) = var(stim_counts_ff, 0) / mu_ff;
        else
            ff(u) = NaN;
        end
    end

    % Save condition-level outputs.
    run_metrics.conditions(c).stim_end_min = stim_end_min;
    run_metrics.conditions(c).stim_end_median = stim_end_med;
    run_metrics.conditions(c).fr_stim = fr_stim;
    run_metrics.conditions(c).dprime = dprime;
    run_metrics.conditions(c).fano_factor = ff;
    run_metrics.conditions(c).fr_stim_trial = fr_stim_trial_mat;
    run_metrics.conditions(c).dprime_stim_rate_trial = dprime_stim_rate_trial_mat;
    run_metrics.conditions(c).prestim_rate_trial = pre_rate_trial_mat;
    run_metrics.conditions(c).ff_spikecount_trial = ff_count_trial_mat;

    % Fill unit-by-condition summary tables.
    fr_by_cond(:, c) = fr_stim;
    dp_by_cond(:, c) = dprime;
    ff_by_cond(:, c) = ff;
end

run_metrics.metric_by_condition.fr_stim = fr_by_cond;
run_metrics.metric_by_condition.dprime = dp_by_cond;
run_metrics.metric_by_condition.fano_factor = ff_by_cond;

% For plotting, keep across-condition summary values for each unit.
run_metrics.max_fr_stim = max(fr_by_cond, [], 2);
run_metrics.min_fr_stim = min(fr_by_cond, [], 2);
run_metrics.max_dprime = max(dp_by_cond, [], 2);
run_metrics.min_dprime = min(dp_by_cond, [], 2);
run_metrics.max_fano_factor = max(ff_by_cond, [], 2);

end


function [condIndex, condDefs, used_fields] = group_trials_by_condition(stiminfo_run)
%% =========================================================================
% group_trials_by_condition
%
% Purpose:
%   Group trials into conditions by comparing all valid trial-wise fields
%   in stiminfo_run, excluding condition_ids.
%
% Inputs:
%   stiminfo_run : struct
%       One run of stimulus information containing trial-wise fields.
%
% Outputs:
%   condIndex : numeric vector
%       Condition index for each trial.
%   condDefs : struct array
%       One struct per condition, containing trial_indices and condition
%       field values.
%   used_fields : cell array of char
%       Trial-wise stiminfo fields used to define conditions.
%
% =========================================================================

% Infer trial count from stiminfo.
nTrial = get_n_trial_from_stiminfo(stiminfo_run);
fn = fieldnames(stiminfo_run);

used_fields = {};

% Select only fields that appear to contain one value per trial.
for i = 1:numel(fn)
    f = fn{i};

    % Explicitly ignore condition_ids.
    if strcmp(f, 'condition_ids')
        continue;
    end

    val = stiminfo_run.(f);

    if is_trialwise_field(val, nTrial)
        used_fields{end+1} = f; %#ok<AGROW>
    end
end

if isempty(used_fields)
    error('No valid trial-wise fields found in stiminfo for condition grouping.');
end

% Build one stable string key per trial.
condKeys = strings(nTrial, 1);

for t = 1:nTrial
    parts = strings(1, numel(used_fields));

    for k = 1:numel(used_fields)
        f = used_fields{k};
        this_val = get_trial_value(stiminfo_run.(f), t, nTrial);
        parts(k) = sprintf('%s=%s', f, serialize_value(this_val));
    end

    condKeys(t) = strjoin(cellstr(parts), ' | ');
end

% unique(...,'stable') preserves first appearance order.
[uniqueKeys, ~, condIndex] = unique(condKeys, 'stable'); %#ok<ASGLU>
nCond = numel(uniqueKeys);

condDefs = repmat(struct(), nCond, 1);

for c = 1:nCond
    idx = find(condIndex == c);
    first_idx = idx(1);

    condDefs(c).trial_indices = idx(:);

    % Store representative field values from the first trial in the group.
    for k = 1:numel(used_fields)
        f = used_fields{k};
        condDefs(c).(f) = get_trial_value(stiminfo_run.(f), first_idx, nTrial);
    end
end

end


function tf = is_trialwise_field(val, nTrial)
%% =========================================================================
% is_trialwise_field
%
% Purpose:
%   Determine whether a stiminfo field appears to store one value per trial.
%
% Inputs:
%   val : any type
%       One field value from stiminfo_run.
%   nTrial : numeric scalar
%       Expected number of trials.
%
% Outputs:
%   tf : logical scalar
%       True if the field is interpreted as trial-wise; false otherwise.
%
% =======================================================================

if iscell(val)
    tf = numel(val) == nTrial;
    return;
end

if isstring(val)
    tf = (numel(val) == nTrial) || (size(val,1) == nTrial);
    return;
end

if ischar(val)
    % A constant char array should not be treated as trial-wise.
    % Only treat it as trial-wise when rows correspond to trials.
    tf = (size(val,1) == nTrial) && (nTrial > 1);
    return;
end

if isnumeric(val) || islogical(val)
    tf = (numel(val) == nTrial) || (size(val,1) == nTrial);
    return;
end

tf = false;
end


function v = get_trial_value(val, idx, nTrial)
%% =========================================================================
% get_trial_value
%
% Purpose:
%   Extract the value of one trial from a trial-wise stiminfo field.
%
% Inputs:
%   val : any supported type
%       Trial-wise field value.
%   idx : numeric scalar
%       Trial index to extract.
%   nTrial : numeric scalar
%       Total number of trials.
%
% Outputs:
%   v : extracted value
%       Value corresponding to the requested trial.
%
% =========================================================================

if iscell(val)
    v = val{idx};
    return;
end

if isstring(val)
    if isvector(val) && numel(val) == nTrial
        v = val(idx);
    elseif size(val,1) == nTrial
        v = val(idx, :);
    else
        error('Unsupported string field format.');
    end
    return;
end

if ischar(val)
    if size(val,1) == nTrial && nTrial > 1
        v = val(idx, :);
    else
        error('Char field is not trial-wise.');
    end
    return;
end

if isnumeric(val) || islogical(val)
    if isvector(val) && numel(val) == nTrial
        v = val(idx);
    elseif size(val,1) == nTrial
        v = val(idx, :);
    else
        error('Numeric/logical field is not trial-wise.');
    end
    return;
end

error('Unsupported field type for trial value extraction.');
end


function s = serialize_value(v)
%% =========================================================================
% serialize_value
%
% Purpose:
%   Convert a trial value into a stable string representation so that it
%   can be used in condition grouping keys.
%
% Inputs:
%   v : any supported type
%       Trial value to serialize.
%
% Outputs:
%   s : string
%       Stable string representation of the input value.
%
% =========================================================================

if isstring(v)
    s = strjoin(cellstr(v(:)), ',');
    s = string(s);
    return;
end

if ischar(v)
    s = string(v);
    return;
end

if isnumeric(v) || islogical(v)
    v = v(:)';
    parts = strings(1, numel(v));

    for i = 1:numel(v)
        if isnan(v(i))
            parts(i) = "NaN";
        elseif isinf(v(i))
            if v(i) > 0
                parts(i) = "Inf";
            else
                parts(i) = "-Inf";
            end
        else
            parts(i) = string(sprintf('%.12g', v(i)));
        end
    end

    s = "[" + strjoin(cellstr(parts), ',') + "]";
    return;
end

if iscell(v)
    parts = strings(1, numel(v));
    for i = 1:numel(v)
        parts(i) = serialize_value(v{i});
    end
    s = "{" + strjoin(cellstr(parts), ',') + "}";
    return;
end

s = string(mat2str(v));
end


function nTrial = get_n_trial_from_stiminfo(stiminfo_run)

%% =========================================================================
% get_n_trial_from_stiminfo
%
% Purpose:
%   Infer the number of trials from the fields stored in a stiminfo struct.
%
% Inputs:
%   stiminfo_run : struct
%       One run of stimulus information.
%
% Outputs:
%   nTrial : numeric scalar
%       Inferred number of trials.
%
% =========================================================================

fn = fieldnames(stiminfo_run);
nCandidates = [];

% Collect candidate trial counts from fields that look trial-wise.
for i = 1:numel(fn)
    val = stiminfo_run.(fn{i});

    if iscell(val)
        if numel(val) > 1
            nCandidates(end+1) = numel(val); %#ok<AGROW>
        end

    elseif isstring(val)
        if isvector(val)
            if numel(val) > 1
                nCandidates(end+1) = numel(val); %#ok<AGROW>
            end
        else
            if size(val,1) > 1
                nCandidates(end+1) = size(val,1); %#ok<AGROW>
            end
        end

    elseif isnumeric(val) || islogical(val)
        if isvector(val)
            if numel(val) > 1
                nCandidates(end+1) = numel(val); %#ok<AGROW>
            end
        else
            if size(val,1) > 1
                nCandidates(end+1) = size(val,1); %#ok<AGROW>
            end
        end
    end
end

if isempty(nCandidates)
    error('Cannot infer trial count from stiminfo.');
end

% Use the mode as a robust estimate of trial count.
nTrial = mode(nCandidates);
end


function stim_end = get_trial_stim_end(tr)
%% =========================================================================
% get_trial_stim_end
%
% Purpose:
%   Extract the stimulus end time from one trial matrix.
%
% Inputs:
%   tr : numeric matrix
%       Trial matrix where column 1 is unit/event ID and column 2 is time.
%
% Outputs:
%   stim_end : numeric scalar
%       Stimulus end time, defined here as the second marker time for
%       unit/event ID 2000.
%
% =========================================================================

stim_end = NaN;

if isempty(tr)
    return;
end

% By convention, event marker ID 2000 is used here, and the second marker
% time is treated as stimulus end.
marker_t = tr(tr(:,1) == 2000, 2);

stim_end = marker_t(2);

end


function plot_condition_summary_maxmin(unit_condition_metrics, ksDir)
%% =========================================================================
% plot_condition_summary_maxmin
%
% Purpose:
%   Plot per-run histograms using, for each unit:
%       - maximum FR across conditions
%       - minimum FR across conditions
%       - maximum d-prime across conditions
%       - minimum d-prime across conditions
%       - maximum Fano factor across conditions
%
% Inputs:
%   unit_condition_metrics : cell array
%       One run_metrics struct per run.
%   ksDir : char
%       Kilosort directory used as save location and figure title.
%
% Outputs:
%   None directly.
%   Saves:
%       - condition_metrics_summary_maxmin.png
%       - condition_metrics_summary_maxmin.fig
%
% =========================================================================

nRun = numel(unit_condition_metrics);

% Histogram bin edges
fr_edges = 0:2:60;
dp_edges = -0.5:0.25:4;
ff_edges = 0:0.25:5;

figW = max(1200, 330 * nRun);
figH = 1350;

hfig = figure('Color', 'w', 'Position', [50 50 figW figH]);

for j = 1:nRun

    run_metrics = unit_condition_metrics{j};

    max_fr_plot = run_metrics.max_fr_stim;
    min_fr_plot = run_metrics.min_fr_stim;
    max_dp_plot = run_metrics.max_dprime;
    min_dp_plot = run_metrics.min_dprime;
    ff_plot     = run_metrics.max_fano_factor;

    % Remove NaN values only for plotting.
    max_fr_plot = max_fr_plot(isfinite(max_fr_plot));
    min_fr_plot = min_fr_plot(isfinite(min_fr_plot));
    max_dp_plot = max_dp_plot(isfinite(max_dp_plot));
    min_dp_plot = min_dp_plot(isfinite(min_dp_plot));
    ff_plot     = ff_plot(isfinite(ff_plot));

    % ------------------------------------------------------------
    % Row 1: maximum firing rate across conditions
    % ------------------------------------------------------------
    subplot(5, nRun, j);
    histogram(max_fr_plot, fr_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.30 0.50 0.90], 'EdgeColor', [0.30 0.50 0.90]);
    title(run_metrics.stim_tag, 'Interpreter', 'none');
    xlabel('max FR across conditions (sp/s)');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(max_fr_plot);
    xlim([fr_edges(1) fr_edges(end)]);

    % ------------------------------------------------------------
    % Row 2: minimum firing rate across conditions
    % ------------------------------------------------------------
    subplot(5, nRun, nRun + j);
    histogram(min_fr_plot, fr_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.55 0.55 0.55], 'EdgeColor', [0.55 0.55 0.55]);
    xlabel('min FR across conditions (sp/s)');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(min_fr_plot);
    xlim([fr_edges(1) fr_edges(end)]);

    % ------------------------------------------------------------
    % Row 3: maximum d-prime across conditions
    % ------------------------------------------------------------
    subplot(5, nRun, 2*nRun + j);
    histogram(max_dp_plot, dp_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.25 0.60 0.90], 'EdgeColor', [0.25 0.60 0.90]);
    xlabel('max d-prime across conditions');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(max_dp_plot);
    xlim([dp_edges(1) dp_edges(end)]);

    % ------------------------------------------------------------
    % Row 4: minimum d-prime across conditions
    % ------------------------------------------------------------
    subplot(5, nRun, 3*nRun + j);
    histogram(min_dp_plot, dp_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.65 0.65 0.65], 'EdgeColor', [0.65 0.65 0.65]);
    xlabel('min d-prime across conditions');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(min_dp_plot);
    xlim([dp_edges(1) dp_edges(end)]);

    % ------------------------------------------------------------
    % Row 5: maximum Fano factor across conditions
    % ------------------------------------------------------------
    subplot(5, nRun, 4*nRun + j);
    histogram(ff_plot, ff_edges, 'Normalization', 'probability', ...
        'FaceColor', [0.20 0.70 0.50], 'EdgeColor', [0.20 0.70 0.50]);
    xlabel('max Fano factor across conditions');
    if j == 1
        ylabel('Proportion of units');
    end
    add_mu_median_text(ff_plot);
    xlim([ff_edges(1) ff_edges(end)]);
end

% Use the ks folder path as the figure title.
sgtitle(strrep(ksDir, '\', '\\'), 'Interpreter', 'none');

saveas(hfig, fullfile(ksDir, 'condition_metrics_summary_maxmin.png'));
savefig(hfig, fullfile(ksDir, 'condition_metrics_summary_maxmin.fig'));
close(hfig);

end



function dp = calc_dprime_from_rates(stim_rates, pre_rates)

%% =========================================================================
% calc_dprime_from_rates
%
% Purpose:
%   Compute d-prime from trial-by-trial stimulus and prestimulus firing
%   rates.
%
% Inputs:
%   stim_rates : numeric vector
%       Trial-wise firing rates during the stimulus window.
%   pre_rates : numeric vector
%       Trial-wise firing rates during the prestimulus window.
%
% Outputs:
%   dp : numeric scalar
%       d-prime value. Returns NaN when the value cannot be computed.
%
% =========================================================================

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

function add_mu_median_text(x)
%% =========================================================================
% add_mu_median_text
%
% Purpose:
%   Add sample size, mean, and median text to the current axis.
%
% Inputs:
%   x : numeric vector
%       Data vector used to compute displayed summary statistics.
%
% Outputs:
%   None directly.
%   Adds text annotation to the current axis.
%
% =========================================================================

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


function stimTag_out = normalize_stimtag_list(stimTag_in)
%% =========================================================================
% normalize_stimtag_list
%
% Purpose:
%   Convert a stimTag list into a column cell array of char row vectors.
%
% Inputs:
%   stimTag_in : cell array or string array
%       Input stimTag list.
%
% Outputs:
%   stimTag_out : cell array
%       Normalized stimTag list as cell array of char.
%
% =========================================================================

if isstring(stimTag_in)
    stimTag_out = cellstr(stimTag_in(:));
elseif iscell(stimTag_in)
    stimTag_out = cell(size(stimTag_in));
    for i = 1:numel(stimTag_in)
        stimTag_out{i} = normalize_one_stimtag(stimTag_in{i});
    end
    stimTag_out = stimTag_out(:);
else
    error('stimTag must be a cell array or string array.');
end
end


function s = normalize_one_stimtag(x)

%% =========================================================================
% normalize_one_stimtag
%
% Purpose:
%   Convert one stim tag into a char vector.
%
% Inputs:
%   x : char or scalar string
%       One stim tag.
%
% Outputs:
%   s : char
%       Normalized stim tag.
%
% =========================================================================

if isstring(x)
    if numel(x) ~= 1
        error('Each stim tag must be scalar.');
    end
    s = char(x);

elseif ischar(x)
    s = x;

else
    error('Each stim tag must be char or scalar string.');
end
end