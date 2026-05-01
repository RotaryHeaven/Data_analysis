%% =========================================================================
% model_data_prepar
%
% Purpose:
%   For selected stimTag runs, prepare model-ready trial data from multiple
%   probe-specific kilosort folders.
%
%   For each selected stim_tag:
%     1) Select units from each probe using one of two user-defined methods.
%     2) Read eight data types from bined_data_allruns.
%     3) Keep only selected units in each probe.
%     4) Merge probes along the unit dimension for each trial.
%     5) Remove NaN-containing trials or neurons using user-defined
%        strategies.
%     6) Convert each data type from unit x trial x T into a 1 x Ntrial
%        struct array with fields:
%            - trialId
%            - T
%            - y   (unit x T)
%     7) Group each processed data type again by condition and store
%        condition-wise trial struct arrays.
%
% Inputs required in each kilosort folder:
%   - unit_run_metrics.mat
%   - unit_condition_metrics.mat
%   - bined_data_allruns.mat
%
% Output saved in the common CatGT folder:
%   - model_data_allruns.mat
%
% Output structure:
%   model_data_allruns{i}
%       .stim_tag
%       .analysis_window
%       .bin_size
%       .fr_threshold
%       .ff_threshold
%       .unit_selection_method
%       .nan_trial_strategy
%       .groupd
%       .probe0_usedunit_ids
%       .probe1_usedunit_ids
%       ...
%       .condition_fields
%       .condition_index_per_trial_full
%       .conditions_full
%       .n_trials_full
%       .raw_count
%       .raw_fr
%       .z_within_trial
%       .z_within_condition
%       .z_across_conditions
%       .demean_count_within_trial
%       .demean_fr_within_trial
%       .demean_pooledsd_within_condition
%
%   For each of the eight main data fields above, an additional field is
%   also stored:
%       .<datafield>_by_condition
%
%   Each .<datafield>_by_condition is a 1 x Ncond struct array. Each
%   condition struct contains:
%       .condition_index
%       .trial_indices_full
%       .trial_ids_present
%       .n_trials_full
%       .n_trials_present
%       .trials
%       plus all condition identity fields from conditions_full except
%       trial_indices
%
%   If nan_trial_strategy == 4, the following additional fields are stored:
%       .raw_count_nanmask
%       .raw_fr_nanmask
%       .z_within_trial_nanmask
%       .z_within_condition_nanmask
%       .z_across_conditions_nanmask
%       .demean_count_within_trial_nanmask
%       .demean_fr_within_trial_nanmask
%       .demean_pooledsd_within_condition_nanmask
%
%   If nan_trial_strategy == 4, each mask field also has:
%       .<maskfield>_by_condition
%
%   If nan_trial_strategy == 6, the following additional field patterns are stored:
%       .<datafield>_groupd
%       .<datafield>_probe0_usedunit_ids
%       .<datafield>_probe1_usedunit_ids
%       ...
%
% Unit selection methods:
%   method 1:
%       keep units with
%           unit_run_metrics.fr_stim    >= fr_threshold
%           unit_run_metrics.fano_factor < ff_threshold
%
%   method 2:
%       keep units with
%           unit_condition_metrics.min_fr_stim      >= fr_threshold
%           unit_condition_metrics.max_fano_factor  < ff_threshold
%
% NaN trial strategies:
%   strategy 1:
%       remove a trial from ALL data types if that trial contains a
%       NaN in ANY one of the merged data types
%
%   strategy 2:
%       remove NaN-containing trials independently for each data type,
%       without forcing the outputs to stay trial-aligned
%
%   strategy 3:
%       do not remove any trials, even if NaN is present
%
%   strategy 4:
%       do not remove any trials; replace NaN with 0 and store a mask
%
%   strategy 5:
%       do not remove any trials; remove any neuron that contains NaN in
%       any one of the merged data types, then update probe-specific
%       used unit IDs and groupd accordingly
%
%   strategy 6:
%       do not remove any trials; for each data type independently,
%       remove neurons that contain NaN in that data type only
%
% Notes:
%   1. All probe folders must belong to the same CatGT folder.
%   2. For a given stim_tag, all probes must have the same analysis_window,
%      bin_size, trial count, bin count, and trial identities.
%   3. For within-trial z-score, NaN from zero variance is common;
%      strategy 4 and 6 are often more suitable if this data type is used.
% =========================================================================

clc; clear;
addpath(genpath(fullfile('.', 'expo_tools')));
addpath(genpath(fullfile('.', 'utils')));

%% ----------------------- User parameters -----------------------

probe_ksDirs = { ...
    'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\RafiL001p0120_g1_imec0\kilosort4_10_dedup_phy', ...
    'I:\np_data\RafiL001p0120_g1\catgt_RafiL001p0120_g1\RafiL001p0120_g1_imec1\kilosort4_2_dedup_phy' ...
    };

stimTag = { ...
    '[RFG_coarse2dg_99_4_150isi]', ...
    '[dir12_gpl_2_200isi_fixedphase]', ...
    '_2[Gpl2_2c_2sz_400_2_200isi]'};

fr_threshold = 0.5;
ff_threshold = 5;

% Unit selection method:
%   1 = use unit_run_metrics
%   2 = use unit_condition_metrics
unit_selection_method = 2;

% 1 = global aligned deletion across all data types
% 2 = each data type removes its own NaN trials independently
% 3 = do not remove any trials
% 4 = do not remove any trials; replace NaN with 0 and store a mask
% 5 = do not remove any trials; remove neurons globally across all data types
% 6 = do not remove any trials; remove neurons independently for each data type
nan_trial_strategy = 4;

%% ----------------------- Validate user parameters -----------------------

if ~iscell(probe_ksDirs) || isempty(probe_ksDirs)
    error('probe_ksDirs must be a non-empty cell array of kilosort folder paths.');
end

for p = 1:numel(probe_ksDirs)
    if ~ischar(probe_ksDirs{p}) && ~isstring(probe_ksDirs{p})
        error('probe_ksDirs{%d} must be a char or string path.', p);
    end
    probe_ksDirs{p} = char(probe_ksDirs{p});
end

if ~isscalar(fr_threshold) || ~isnumeric(fr_threshold) || ~isfinite(fr_threshold)
    error('fr_threshold must be a finite numeric scalar.');
end

if ~isscalar(ff_threshold) || ~isnumeric(ff_threshold) || ~isfinite(ff_threshold)
    error('ff_threshold must be a finite numeric scalar.');
end

if ~(isequal(unit_selection_method, 1) || isequal(unit_selection_method, 2))
    error('unit_selection_method must be 1 or 2.');
end

if ~(isequal(nan_trial_strategy, 1) || isequal(nan_trial_strategy, 2) || ...
        isequal(nan_trial_strategy, 3) || isequal(nan_trial_strategy, 4) || ...
        isequal(nan_trial_strategy, 5) || isequal(nan_trial_strategy, 6))
    error('nan_trial_strategy must be 1, 2, 3, 4, 5, or 6.');
end

%% ----------------------- Determine common CatGT folder -----------------------

catgt_folder = get_common_catgt_folder(probe_ksDirs);
fprintf('Common catgt_folder: %s\n', catgt_folder);

%% ----------------------- Process all selected stim tags -----------------------

model_data_allruns = cell(numel(stimTag), 1);

for s = 1:numel(stimTag)

    this_stim_tag = stimTag{s};

    fprintf('\n============================================================\n');
    fprintf('Processing stim_tag: %s\n', this_stim_tag);
    fprintf('============================================================\n');

    probe_data = cell(numel(probe_ksDirs), 1);
    groupd = zeros(1, numel(probe_ksDirs));

    for p = 1:numel(probe_ksDirs)

        ksDir = probe_ksDirs{p};
        fprintf('  Probe %d ksDir: %s\n', p-1, ksDir);

        probe_data{p} = process_one_probe_one_run( ...
            ksDir, this_stim_tag, fr_threshold, ff_threshold, unit_selection_method);

        groupd(p) = numel(probe_data{p}.used_unit_ids);
    end

    ref = probe_data{1};
    for p = 2:numel(probe_ksDirs)
        validate_probe_alignment(ref, probe_data{p}, this_stim_tag, p-1);
    end

    merged = merge_probes_for_one_run(probe_data);

    model_data_allruns{s} = build_model_output_for_one_run( ...
        merged, probe_data, this_stim_tag, fr_threshold, ff_threshold, ...
        unit_selection_method, nan_trial_strategy, groupd);

end

%% ----------------------- Save output -----------------------

save(fullfile(catgt_folder, 'model_data_allruns.mat'), 'model_data_allruns');

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(catgt_folder, 'model_data_allruns.mat'));
fprintf('\nDone.\n');

%% ======================= Local functions =======================

function probe_out = process_one_probe_one_run(ksDir, stim_tag, fr_threshold, ff_threshold, unit_selection_method)

if ~isfolder(ksDir)
    error('kilosort folder does not exist: %s', ksDir);
end

run_file = fullfile(ksDir, 'unit_run_metrics.mat');
cond_file = fullfile(ksDir, 'unit_condition_metrics.mat');
bined_file = fullfile(ksDir, 'bined_data_allruns.mat');

if ~isfile(run_file)
    error('Missing file: %s', run_file);
end
if ~isfile(cond_file)
    error('Missing file: %s', cond_file);
end
if ~isfile(bined_file)
    error('Missing file: %s', bined_file);
end

Srun = load(run_file, 'unit_run_metrics');
Scond = load(cond_file, 'unit_condition_metrics');
Sbined = load(bined_file, 'bined_data_allruns');

if ~isfield(Srun, 'unit_run_metrics')
    error('unit_run_metrics not found in %s', run_file);
end
if ~isfield(Scond, 'unit_condition_metrics')
    error('unit_condition_metrics not found in %s', cond_file);
end
if ~isfield(Sbined, 'bined_data_allruns')
    error('bined_data_allruns not found in %s', bined_file);
end

unit_run_metrics = Srun.unit_run_metrics;
unit_condition_metrics = Scond.unit_condition_metrics;
bined_data_allruns = Sbined.bined_data_allruns;

run_idx_in_run_metrics = find_run_index_by_stim_tag(unit_run_metrics, stim_tag);
run_idx_in_cond_metrics = find_run_index_by_stim_tag(unit_condition_metrics, stim_tag);
run_idx_in_bined = find_run_index_by_stim_tag(bined_data_allruns, stim_tag);

run_metrics = unit_run_metrics{run_idx_in_run_metrics};
cond_metrics = unit_condition_metrics{run_idx_in_cond_metrics};
bined_entry = bined_data_allruns{run_idx_in_bined};

[used_unit_ids, ~] = select_units( ...
    run_metrics, cond_metrics, fr_threshold, ff_threshold, unit_selection_method);

if ~isfield(bined_entry, 'unit_ids')
    error('unit_ids missing in bined_data_allruns entry for stim_tag %s', stim_tag);
end

bined_unit_ids = bined_entry.unit_ids(:);
[tf, idx_in_bined] = ismember(used_unit_ids, bined_unit_ids);

if ~all(tf)
    missing_ids = used_unit_ids(~tf);
    error('Selected unit_ids not found in bined_data_allruns for stim_tag %s: %s', ...
        stim_tag, mat2str(missing_ids(:)'));
end

data_fields = get_data_field_list();

probe_out = struct();
probe_out.stim_tag = stim_tag;
probe_out.analysis_window = bined_entry.analysis_window;
probe_out.bin_size = bined_entry.bin_size;
probe_out.bin_edges = bined_entry.bin_edges;
probe_out.bin_centers = bined_entry.bin_centers;

probe_out.unit_ids_all = bined_unit_ids;
probe_out.used_unit_ids = used_unit_ids(:);

probe_out.condition_fields = bined_entry.condition_fields;
probe_out.condition_index_per_trial = bined_entry.condition_index_per_trial(:);
probe_out.conditions = bined_entry.conditions;

for k = 1:numel(data_fields)
    f = data_fields{k};
    if ~isfield(bined_entry, f)
        error('Field %s missing in bined_data_allruns entry for stim_tag %s', f, stim_tag);
    end

    X = bined_entry.(f);

    if size(X, 1) ~= numel(bined_unit_ids)
        error('Field %s has unit dimension mismatch in ksDir %s', f, ksDir);
    end

    probe_out.(f) = X(idx_in_bined, :, :);
end

end


function [used_unit_ids, keep_idx] = select_units(run_metrics, cond_metrics, fr_threshold, ff_threshold, unit_selection_method)

switch unit_selection_method
    case 1
        if ~isfield(run_metrics, 'unit_ids') || ~isfield(run_metrics, 'fr_stim') || ~isfield(run_metrics, 'fano_factor')
            error('unit_run_metrics entry is missing required fields for method 1.');
        end

        unit_ids = run_metrics.unit_ids(:);
        fr_metric = run_metrics.fr_stim(:);
        ff_metric = run_metrics.fano_factor(:);

        keep_idx = isfinite(fr_metric) & isfinite(ff_metric) & ...
            (fr_metric >= fr_threshold) & (ff_metric < ff_threshold);

        used_unit_ids = unit_ids(keep_idx);

    case 2
        if ~isfield(cond_metrics, 'unit_ids') || ~isfield(cond_metrics, 'min_fr_stim') || ~isfield(cond_metrics, 'max_fano_factor')
            error('unit_condition_metrics entry is missing required fields for method 2.');
        end

        unit_ids = cond_metrics.unit_ids(:);
        fr_metric = cond_metrics.min_fr_stim(:);
        ff_metric = cond_metrics.max_fano_factor(:);

        keep_idx = isfinite(fr_metric) & isfinite(ff_metric) & ...
            (fr_metric >= fr_threshold) & (ff_metric < ff_threshold);

        used_unit_ids = unit_ids(keep_idx);

    otherwise
        error('Unknown unit_selection_method.');
end

end


function validate_probe_alignment(ref, cur, stim_tag, probe_index)

if ~isequal(ref.analysis_window, cur.analysis_window)
    error('analysis_window mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

if ~isequal(ref.bin_size, cur.bin_size)
    error('bin_size mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

if ~isequal(ref.bin_edges, cur.bin_edges)
    error('bin_edges mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

if ~isequal(ref.bin_centers, cur.bin_centers)
    error('bin_centers mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

if ~isequal(ref.condition_index_per_trial, cur.condition_index_per_trial)
    error('condition_index_per_trial mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

if numel(ref.conditions) ~= numel(cur.conditions)
    error('Condition count mismatch across probes for stim_tag %s (probe %d).', ...
        stim_tag, probe_index);
end

data_fields = get_data_field_list();

for k = 1:numel(data_fields)
    f = data_fields{k};

    Xref = ref.(f);
    Xcur = cur.(f);

    if size(Xref, 2) ~= size(Xcur, 2)
        error('Trial count mismatch in %s across probes for stim_tag %s (probe %d).', ...
            f, stim_tag, probe_index);
    end

    if size(Xref, 3) ~= size(Xcur, 3)
        error('Bin count mismatch in %s across probes for stim_tag %s (probe %d).', ...
            f, stim_tag, probe_index);
    end
end

end


function merged = merge_probes_for_one_run(probe_data)

data_fields = get_data_field_list();

merged = struct();
merged.stim_tag = probe_data{1}.stim_tag;
merged.analysis_window = probe_data{1}.analysis_window;
merged.bin_size = probe_data{1}.bin_size;
merged.bin_edges = probe_data{1}.bin_edges;
merged.bin_centers = probe_data{1}.bin_centers;

merged.condition_fields = probe_data{1}.condition_fields;
merged.condition_index_per_trial_full = probe_data{1}.condition_index_per_trial;
merged.conditions_full = probe_data{1}.conditions;

for k = 1:numel(data_fields)
    f = data_fields{k};

    Xcat = probe_data{1}.(f);
    for p = 2:numel(probe_data)
        Xcat = cat(1, Xcat, probe_data{p}.(f));
    end
    merged.(f) = Xcat;
end

end


function out = build_model_output_for_one_run(merged, probe_data, stim_tag, fr_threshold, ff_threshold, unit_selection_method, nan_trial_strategy, groupd)

data_fields = get_data_field_list();

out = struct();
out.stim_tag = stim_tag;
out.analysis_window = merged.analysis_window;
out.bin_size = merged.bin_size;

out.fr_threshold = fr_threshold;
out.ff_threshold = ff_threshold;
out.unit_selection_method = unit_selection_method;
out.nan_trial_strategy = nan_trial_strategy;

if nan_trial_strategy ~= 6
    out.groupd = groupd(:)';

    for p = 1:numel(probe_data)
        field_name = sprintf('probe%d_usedunit_ids', p-1);
        out.(field_name) = probe_data{p}.used_unit_ids(:);
    end
end

out.condition_fields = merged.condition_fields;
out.condition_index_per_trial_full = merged.condition_index_per_trial_full;
out.conditions_full = merged.conditions_full;
out.n_trials_full = numel(merged.condition_index_per_trial_full);

switch nan_trial_strategy
    case 1
        bad_trial = false(1, out.n_trials_full);

        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_trial = bad_trial | get_bad_trial_mask(merged.(f));
        end

        keep_trial_ids = find(~bad_trial);
        out.kept_trial_ids_global = keep_trial_ids(:)';

        for k = 1:numel(data_fields)
            f = data_fields{k};
            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 2
        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_trial = get_bad_trial_mask(merged.(f));
            keep_trial_ids = find(~bad_trial);

            keep_field_name = sprintf('%s_kept_trial_ids', f);
            out.(keep_field_name) = keep_trial_ids(:)';

            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 3
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;

        for k = 1:numel(data_fields)
            f = data_fields{k};
            out.(f) = build_trial_struct_array(merged.(f), keep_trial_ids);
        end

    case 4
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;

        for k = 1:numel(data_fields)
            f = data_fields{k};

            X = merged.(f);
            nanmask = isnan(X);

            X_filled = X;
            X_filled(nanmask) = 0;

            out.(f) = build_trial_struct_array(X_filled, keep_trial_ids);

            mask_field_name = sprintf('%s_nanmask', f);
            out.(mask_field_name) = build_trial_struct_array(nanmask, keep_trial_ids);
        end

    case 5
        bad_neuron = false(sum(groupd), 1);

        for k = 1:numel(data_fields)
            f = data_fields{k};
            bad_neuron = bad_neuron | get_bad_neuron_mask(merged.(f));
        end

        keep_neuron = ~bad_neuron;

        if ~any(keep_neuron)
            error('After nan_trial_strategy = 5, no neurons remain for stim_tag %s.', stim_tag);
        end

        [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal( ...
            probe_data, groupd, keep_neuron);

        out.groupd = new_groupd(:)';

        for p = 1:numel(new_probe_unit_ids)
            field_name = sprintf('probe%d_usedunit_ids', p-1);
            out.(field_name) = new_probe_unit_ids{p}(:);
        end

        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;
        out.kept_neuron_global = find(keep_neuron(:))';

        for k = 1:numel(data_fields)
            f = data_fields{k};
            X = merged.(f);
            X = X(keep_neuron, :, :);
            out.(f) = build_trial_struct_array(X, keep_trial_ids);
        end

    case 6
        keep_trial_ids = 1:out.n_trials_full;
        out.kept_trial_ids_global = keep_trial_ids;

        for k = 1:numel(data_fields)
            f = data_fields{k};

            X = merged.(f);
            bad_neuron = get_bad_neuron_mask(X);
            keep_neuron = ~bad_neuron;

            [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal( ...
                probe_data, groupd, keep_neuron);

            groupd_field = sprintf('%s_groupd', f);
            out.(groupd_field) = new_groupd(:)';

            for p = 1:numel(new_probe_unit_ids)
                unit_field_name = sprintf('%s_probe%d_usedunit_ids', f, p-1);
                out.(unit_field_name) = new_probe_unit_ids{p}(:);
            end

            kept_neuron_field = sprintf('%s_kept_neuron_global', f);
            out.(kept_neuron_field) = find(keep_neuron(:))';

            X = X(keep_neuron, :, :);
            out.(f) = build_trial_struct_array(X, keep_trial_ids);
        end

    otherwise
        error('Unknown nan_trial_strategy.');
end

% -------------------------------------------------------------------------
% Group all processed main data fields by condition
% -------------------------------------------------------------------------
out = add_condition_groupings_to_output(out, data_fields);

% -------------------------------------------------------------------------
% If strategy 4 is used, group all mask fields by condition as well
% -------------------------------------------------------------------------
if nan_trial_strategy == 4
    mask_fields = get_mask_field_list(data_fields);
    out = add_condition_groupings_to_output(out, mask_fields);
end

end


function out = add_condition_groupings_to_output(out, fields_to_group)
%% =========================================================================
% add_condition_groupings_to_output
%
% Purpose:
%   For each specified field already stored in "out", add a new field
%   named "<field>_by_condition", which reorganizes that trial-struct data
%   into condition-wise groups.
% =========================================================================

for k = 1:numel(fields_to_group)
    f = fields_to_group{k};

    if ~isfield(out, f)
        error('Field %s is missing when trying to add condition groupings.', f);
    end

    by_field_name = sprintf('%s_by_condition', f);
    out.(by_field_name) = build_by_condition_struct( ...
        out.(f), out.conditions_full, out.condition_index_per_trial_full);
end

end


function cond_struct = build_by_condition_struct(trial_struct_array, conditions_full, condition_index_per_trial_full)
%% =========================================================================
% build_by_condition_struct
%
% Purpose:
%   Reorganize one processed data field from a 1 x Ntrial struct array into
%   a 1 x Ncond struct array.
%
% Inputs:
%   trial_struct_array : struct array
%       Trial-wise data with fields:
%           - trialId
%           - T
%           - y
%   conditions_full : struct array
%       Full condition definitions for the run.
%   condition_index_per_trial_full : numeric vector
%       Full-length condition label for each original trial index.
%
% Output:
%   cond_struct : struct array
%       Condition-wise grouped data.
% =========================================================================

nCond = numel(conditions_full);
cond_struct = repmat(struct(), 1, nCond);

% Extract all currently present trial IDs from this processed field
present_trial_ids = zeros(1, numel(trial_struct_array));
for i = 1:numel(trial_struct_array)
    present_trial_ids(i) = trial_struct_array(i).trialId;
end

for c = 1:nCond

    % Full trial indices belonging to this condition before any strategy-specific filtering
    if isfield(conditions_full(c), 'trial_indices')
        trial_indices_full = conditions_full(c).trial_indices(:)';
    else
        trial_indices_full = find(condition_index_per_trial_full == c);
    end

    % Keep only trial structs whose trialId maps to this condition
    keep_mask = false(1, numel(trial_struct_array));

    for i = 1:numel(trial_struct_array)
        tid = trial_struct_array(i).trialId;

        if tid < 1 || tid > numel(condition_index_per_trial_full)
            error('trialId %d is out of range for condition_index_per_trial_full.', tid);
        end

        keep_mask(i) = (condition_index_per_trial_full(tid) == c);
    end

    trial_ids_present = present_trial_ids(keep_mask);
    trials_present = trial_struct_array(keep_mask);

    % Basic summary fields
    cond_struct(c).condition_index = c;
    cond_struct(c).trial_indices_full = trial_indices_full(:)';
    cond_struct(c).trial_ids_present = trial_ids_present(:)';
    cond_struct(c).n_trials_full = numel(trial_indices_full);
    cond_struct(c).n_trials_present = numel(trial_ids_present);
    cond_struct(c).trials = trials_present;

    % Copy all condition identity fields except trial_indices
    fn = fieldnames(conditions_full(c));
    for j = 1:numel(fn)
        this_field = fn{j};
        if strcmp(this_field, 'trial_indices')
            continue;
        end
        cond_struct(c).(this_field) = conditions_full(c).(this_field);
    end
end

end


function mask_fields = get_mask_field_list(data_fields)
%% =========================================================================
% get_mask_field_list
%
% Purpose:
%   Build the list of corresponding mask field names for strategy 4.
% =========================================================================

mask_fields = cell(size(data_fields));
for k = 1:numel(data_fields)
    mask_fields{k} = sprintf('%s_nanmask', data_fields{k});
end

end


function bad_trial = get_bad_trial_mask(X)
bad_trial = squeeze(any(any(isnan(X), 1), 3));
bad_trial = reshape(bad_trial, 1, []);
end


function S = build_trial_struct_array(X, keep_trial_ids)

nKeep = numel(keep_trial_ids);
nUnit = size(X, 1);
nBin = size(X, 3);

S = repmat(struct('trialId', [], 'T', [], 'y', []), 1, nKeep);

for i = 1:nKeep
    tr = keep_trial_ids(i);

    S(i).trialId = tr;
    S(i).T = nBin;
    S(i).y = reshape(X(:, tr, :), nUnit, nBin);
end

end


function idx = find_run_index_by_stim_tag(cell_of_structs, stim_tag)

all_tags = cell(numel(cell_of_structs), 1);

for i = 1:numel(cell_of_structs)
    if ~isfield(cell_of_structs{i}, 'stim_tag')
        error('Entry %d is missing stim_tag.', i);
    end
    all_tags{i} = cell_of_structs{i}.stim_tag;
end

idx = find(strcmp(all_tags, stim_tag));

if isempty(idx)
    error('Requested stim_tag not found: %s', stim_tag);
end

if numel(idx) > 1
    error('Duplicate stim_tag found: %s', stim_tag);
end

end


function fields = get_data_field_list()

fields = { ...
    'raw_count', ...
    'raw_fr', ...
    'z_within_trial', ...
    'z_within_condition', ...
    'z_across_conditions', ...
    'demean_count_within_trial', ...
    'demean_fr_within_trial', ...
    'demean_pooledsd_within_condition'};

end


function catgt_folder = get_common_catgt_folder(probe_ksDirs)

catgt_list = cell(numel(probe_ksDirs), 1);

for i = 1:numel(probe_ksDirs)
    ksDir = probe_ksDirs{i};

    probe_folder = fileparts(ksDir);
    catgt_folder_i = fileparts(probe_folder);

    catgt_list{i} = catgt_folder_i;
end

catgt_folder = catgt_list{1};

for i = 2:numel(catgt_list)
    if ~strcmp(catgt_folder, catgt_list{i})
        error(['All probe kilosort folders must belong to the same CatGT folder.\n' ...
            'Got:\n%s\n%s'], catgt_folder, catgt_list{i});
    end
end

end


function bad_neuron = get_bad_neuron_mask(X)

bad_neuron = squeeze(any(any(isnan(X), 2), 3));
bad_neuron = bad_neuron(:);

end


function [new_probe_unit_ids, new_groupd] = update_probe_unit_ids_after_neuron_removal(probe_data, old_groupd, keep_neuron)

new_probe_unit_ids = cell(numel(probe_data), 1);
new_groupd = zeros(1, numel(probe_data));

row_start = 1;

for p = 1:numel(probe_data)
    row_end = row_start + old_groupd(p) - 1;

    this_keep = keep_neuron(row_start:row_end);
    this_ids = probe_data{p}.used_unit_ids(:);

    new_probe_unit_ids{p} = this_ids(this_keep);
    new_groupd(p) = numel(new_probe_unit_ids{p});

    row_start = row_end + 1;
end

end