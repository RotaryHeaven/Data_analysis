%% subspace_overlap_sigtest.m
% Matched-control significance tests for DLAG subspace overlap.
%
% Run this script from one catgt_* session folder after
% subspace_similarity_dlag.m has saved:
%   FA_Dlag_<data_content>/mat_results/runXXX/
%       subspace_similarity_<latent-selection-tag>.mat
%
% For every neural group and for both comparisons
%   1) Across vs Within
%   2) Feedforward vs Feedback
% this script uses one matched random draw to calculate all of:
%   - the full ordered principal-angle spectrum;
%   - A captures B directional similarity;
%   - B captures A directional similarity.
%
% The two available controls are:
%   'uniform'
%       Uniform (Haar) random subspaces in the p-dimensional neural space.
%   'covariance_weighted'
%       Random vectors drawn from N(0,Sigma), where Sigma is calculated
%       separately for each neural group from every trial and time bin in
%       model_data_allruns{run}.(data_content). No DLAG refitting is done.
%
% The raw p-by-(k+m) random matrix is split into its k- and m-column
% blocks first. Each block and the joint matrix must be full rank. The two
% blocks are then orthonormalized separately; the two resulting subspaces
% are NOT forced to be mutually orthogonal.
%
% Outputs:
%   subspace_overlap_sigtest_<control>_<latent-selection-tag>.mat
%   ..._principal_angles.csv
%   ..._principal_angle_fractions.csv
%   ..._similarity.csv

clc;
clear;

%% ========================== USER SETTINGS ==============================

data_content = 'raw_count';
% Common options:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial,
% demean_fr_within_trial, demean_pooledsd_within_condition

runIdx = 1;

% Used only by covariance_weighted control to select the same run from
% model_data_allruns. This is the stimulus tag, not runIdx above.
stim_tag = '[Gpl2_2c_2sz_400_2_200isi]';
dat_file = fullfile('.', 'model_data_allruns');

% Display/file labels only. Order must match SubspaceSim.group order.
group_names = {'V1', 'MT'};

% These settings must match subspace_similarity_dlag.m.
use_dsl_filter = false;
dsl_field = 'logical';             % usually 'rawlogical' or 'logical'

use_svexp_filter = false;
svexp_field = 'logical';           % 'rawlogical' or 'logical'
shared_varexp_threshold = 0.95;

% Control choice:
%   'uniform'
%   'covariance_weighted'
control_method = 'covariance_weighted';

num_null = 10000;
alpha = 0.05;
random_seed = 1;
max_draw_attempts = 1000;
progress_every = 1000;

save_results = true;
save_tables = true;

%% ======================== VALIDATE SETTINGS ============================

group_names = normalizeGroupNamesLocal(group_names);
[group_display_names, group_file_tags] = buildGroupLabelsLocal(group_names);
num_groups = numel(group_names);

control_method = normalizeControlMethodLocal(control_method);

validateattributes(runIdx, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'runIdx');
validateattributes(num_null, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'num_null');
validateattributes(alpha, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1}, mfilename, 'alpha');
validateattributes(random_seed, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'}, mfilename, 'random_seed');
validateattributes(max_draw_attempts, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'max_draw_attempts');
validateattributes(progress_every, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'}, ...
    mfilename, 'progress_every');

latent_selection_tag = makeLatentSelectionTagLocal( ...
    use_dsl_filter, dsl_field, ...
    use_svexp_filter, svexp_field, shared_varexp_threshold);

base_dir = fullfile('.', ['FA_Dlag_', data_content]);
result_dir = fullfile(base_dir, 'mat_results', sprintf('run%03d', runIdx));
input_file = fullfile(result_dir, ...
    sprintf('subspace_similarity_%s.mat', latent_selection_tag));

if ~isfile(input_file)
    error(['Input file not found: %s\nRun subspace_similarity_dlag.m with ', ...
           'the same data/filter settings first.'], input_file);
end

fprintf('Input subspace file : %s\n', input_file);
fprintf('Control method      : %s\n', control_method);
fprintf('Number of null draws: %d\n', num_null);
fprintf('Alpha               : %.6g\n', alpha);
fprintf('Group labels        : %s\n', strjoin(group_display_names, ' | '));

%% ======================= LOAD OBSERVED RESULTS =========================

S = load(input_file, 'SubspaceSim');
if ~isfield(S, 'SubspaceSim')
    error('%s does not contain SubspaceSim.', input_file);
end
SubspaceSim = S.SubspaceSim;
validateSubspaceSimLocal(SubspaceSim, num_groups);

% p is inferred from the saved bases. Covariance-weighted control also
% cross-checks it against model_data_allruns yDims.
group_p = nan(1, num_groups);
for g = 1:num_groups
    group_p(g) = inferGroupObservationDimLocal(SubspaceSim, g);
end

covariance_info = cell(1, num_groups);
sampling_factor = cell(1, num_groups);

if strcmp(control_method, 'covariance_weighted')
    [covariance_info, sampling_factor] = loadGroupCovariancesLocal( ...
        dat_file, stim_tag, data_content, group_p, num_groups);
else
    for g = 1:num_groups
        covariance_info{g} = struct( ...
            'used', false, ...
            'matrix', [], ...
            'num_samples', NaN, ...
            'source_file', '', ...
            'stim_tag', '', ...
            'data_content', '');
        sampling_factor{g} = [];
    end
end

%% ======================= GENERATE MATCHED NULLS ========================

rng(random_seed, 'twister');

pair_names = {'across_vs_within', 'feedforward_vs_feedback'};

SubspaceOverlapSigTest = struct();
SubspaceOverlapSigTest.meta = struct();
SubspaceOverlapSigTest.meta.program = mfilename;
SubspaceOverlapSigTest.meta.created = datestr(now, 30);
SubspaceOverlapSigTest.meta.model_mode = 'all_condition_model';
SubspaceOverlapSigTest.meta.data_content = data_content;
SubspaceOverlapSigTest.meta.runIdx = runIdx;
SubspaceOverlapSigTest.meta.stim_tag = stim_tag;
SubspaceOverlapSigTest.meta.control_method = control_method;
SubspaceOverlapSigTest.meta.num_null = num_null;
SubspaceOverlapSigTest.meta.alpha = alpha;
SubspaceOverlapSigTest.meta.random_seed = random_seed;
SubspaceOverlapSigTest.meta.latent_selection_tag = latent_selection_tag;
SubspaceOverlapSigTest.meta.use_dsl_filter = logical(use_dsl_filter);
SubspaceOverlapSigTest.meta.dsl_field = dsl_field;
SubspaceOverlapSigTest.meta.use_svexp_filter = logical(use_svexp_filter);
SubspaceOverlapSigTest.meta.svexp_field = svexp_field;
SubspaceOverlapSigTest.meta.shared_varexp_threshold = shared_varexp_threshold;
SubspaceOverlapSigTest.meta.group_names = group_names;
SubspaceOverlapSigTest.meta.group_display_names = group_display_names;
SubspaceOverlapSigTest.meta.group_file_tags = group_file_tags;
SubspaceOverlapSigTest.meta.input_file = input_file;
SubspaceOverlapSigTest.meta.model_data_file = dat_file;
SubspaceOverlapSigTest.meta.covariance_source = ...
    'all trials and all time bins in model_data_allruns.(data_content), separately by neural group';
SubspaceOverlapSigTest.meta.null_pairing_note = ...
    ['Within each pair and null iteration, principal angles and both ', ...
     'directional similarities come from the same random subspace draw.'];

principal_rows = {};
fraction_rows = {};
similarity_rows = {};

for g = 1:num_groups
    fprintf('\n============================================================\n');
    fprintf('%s (p = %d)\n', group_display_names{g}, group_p(g));

    group_out = struct();
    group_out.name = group_display_names{g};
    group_out.file_tag = group_file_tags{g};
    group_out.p = group_p(g);
    group_out.covariance = covariance_info{g};
    group_out.pairNames = pair_names;
    group_out.pair = cell(1, numel(pair_names));

    for pair_idx = 1:numel(pair_names)
        pair_name = pair_names{pair_idx};
        observed_pair = getPairByNameLocal(SubspaceSim, g, pair_name);

        pair_out = initializePairOutputLocal(observed_pair, pair_name);

        if ~strcmp(pair_out.status, 'ok')
            fprintf('  %-28s skipped: %s\n', pair_name, pair_out.warning);
            group_out.pair{pair_idx} = pair_out;
            continue;
        end

        QA_obs = observed_pair.basisA;
        QB_obs = observed_pair.basisB;
        p = size(QA_obs, 1);
        k = size(QA_obs, 2);
        m = size(QB_obs, 2);
        r = min(k, m);

        if size(QB_obs, 1) ~= p
            error('%s %s: basisA and basisB have different row counts.', ...
                group_display_names{g}, pair_name);
        end
        if p ~= group_p(g)
            error('%s %s: saved basis has p=%d, expected p=%d.', ...
                group_display_names{g}, pair_name, p, group_p(g));
        end
        if k + m >= p
            error('%s %s violates the required k+m<p condition (%d+%d >= %d).', ...
                group_display_names{g}, pair_name, k, m, p);
        end
        if rank(QA_obs) ~= k || rank(QB_obs) ~= m
            error('%s %s: a saved observed basis is not full column rank.', ...
                group_display_names{g}, pair_name);
        end

        field_ab = makeCaptureFieldNameLocal( ...
            observed_pair.labelA, observed_pair.labelB);
        field_ba = makeCaptureFieldNameLocal( ...
            observed_pair.labelB, observed_pair.labelA);

        observed_angles = reshape(double( ...
            observed_pair.principal.angle_deg), 1, []);
        if numel(observed_angles) ~= r || any(~isfinite(observed_angles))
            error('%s %s: observed angle spectrum must contain %d finite values.', ...
                group_display_names{g}, pair_name, r);
        end

        observed_sim_ab = getObservedSimilarityLocal(observed_pair, field_ab);
        observed_sim_ba = getObservedSimilarityLocal(observed_pair, field_ba);

        null_angles = nan(num_null, r);
        null_sim_ab = nan(num_null, 1);
        null_sim_ba = nan(num_null, 1);

        fprintf('  %-28s k=%d, m=%d, r=%d\n', pair_name, k, m, r);

        for b = 1:num_null
            [QA_null, QB_null] = drawMatchedSubspacesLocal( ...
                p, k, m, control_method, sampling_factor{g}, ...
                max_draw_attempts);

            [this_angles, this_sim_ab, this_sim_ba] = ...
                calculateAllOverlapMetricsLocal(QA_null, QB_null);

            null_angles(b, :) = this_angles;
            null_sim_ab(b) = this_sim_ab;
            null_sim_ba(b) = this_sim_ba;

            if progress_every > 0 && ...
                    (mod(b, progress_every) == 0 || b == num_null)
                fprintf('    null draw %d/%d\n', b, num_null);
            end
        end

        principal_test = principalAngleTestLocal( ...
            observed_angles, null_angles, alpha);
        similarity_test_ab = scalarControlTestLocal( ...
            observed_sim_ab, null_sim_ab, alpha);
        similarity_test_ba = scalarControlTestLocal( ...
            observed_sim_ba, null_sim_ba, alpha);

        pair_out.p = p;
        pair_out.k = k;
        pair_out.m = m;
        pair_out.r = r;
        pair_out.capture_field_A_captures_B = field_ab;
        pair_out.capture_field_B_captures_A = field_ba;

        pair_out.observed = struct();
        pair_out.observed.principal_angle_deg = observed_angles;
        pair_out.observed.similarity = struct();
        pair_out.observed.similarity.(field_ab) = observed_sim_ab;
        pair_out.observed.similarity.(field_ba) = observed_sim_ba;

        pair_out.null = struct();
        pair_out.null.draw_id = (1:num_null)';
        pair_out.null.principal_angle_deg = null_angles;
        pair_out.null.similarity = struct();
        pair_out.null.similarity.(field_ab) = null_sim_ab;
        pair_out.null.similarity.(field_ba) = null_sim_ba;

        pair_out.principal = principal_test;
        pair_out.similarity = struct();
        pair_out.similarity.(field_ab) = similarity_test_ab;
        pair_out.similarity.(field_ba) = similarity_test_ba;

        for j = 1:r
            principal_rows(end + 1, :) = { ... %#ok<SAGROW>
                g, group_display_names{g}, pair_name, ...
                sprintf('%s vs %s', pair_out.labelA, pair_out.labelB), ...
                k, m, r, j, observed_angles(j), ...
                principal_test.null_median_deg(j), ...
                principal_test.p_small(j), ...
                principal_test.p_large(j), ...
                principal_test.p_two(j), ...
                principal_test.direction{j}, ...
                principal_test.significant_aligned(j), ...
                principal_test.significant_orthogonal(j)};
        end

        fraction_rows(end + 1, :) = { ... %#ok<SAGROW>
            g, group_display_names{g}, pair_name, ...
            sprintf('%s vs %s', pair_out.labelA, pair_out.labelB), ...
            k, m, r, principal_test.num_aligned, ...
            principal_test.num_orthogonal, ...
            principal_test.fraction_aligned, ...
            principal_test.fraction_orthogonal, ...
            median(principal_test.null_fraction_aligned), ...
            median(principal_test.null_fraction_orthogonal)};

        similarity_rows(end + 1, :) = makeSimilarityTableRowLocal( ... %#ok<SAGROW>
            g, group_display_names{g}, pair_name, ...
            sprintf('%s captures %s', pair_out.labelA, pair_out.labelB), ...
            k, m, similarity_test_ab);
        similarity_rows(end + 1, :) = makeSimilarityTableRowLocal( ... %#ok<SAGROW>
            g, group_display_names{g}, pair_name, ...
            sprintf('%s captures %s', pair_out.labelB, pair_out.labelA), ...
            k, m, similarity_test_ba);

        group_out.pair{pair_idx} = pair_out;
    end

    SubspaceOverlapSigTest.group(g) = group_out;
end

%% =========================== TABLES / SAVE =============================

principal_variable_names = { ...
    'GroupIndex', 'GroupName', 'PairName', 'Comparison', ...
    'DimA', 'DimB', 'NumAngles', 'AngleIndex', ...
    'ObservedAngleDeg', 'NullMedianDeg', ...
    'PSmall', 'PLarge', 'PTwo', 'Direction', ...
    'SignificantAligned', 'SignificantOrthogonal'};

fraction_variable_names = { ...
    'GroupIndex', 'GroupName', 'PairName', 'Comparison', ...
    'DimA', 'DimB', 'NumAngles', ...
    'NumSignificantAligned', 'NumSignificantOrthogonal', ...
    'ObservedFractionAligned', 'ObservedFractionOrthogonal', ...
    'NullMedianFractionAligned', 'NullMedianFractionOrthogonal'};

similarity_variable_names = { ...
    'GroupIndex', 'GroupName', 'PairName', 'DirectionName', ...
    'DimA', 'DimB', 'ObservedSimilarity', 'NullMedian', 'Effect', ...
    'PHigh', 'PLow', 'PTwo', 'PDirectional', 'Direction', ...
    'SignificantDirectional', 'Star'};

PrincipalAngleTable = cellRowsToTableLocal( ...
    principal_rows, principal_variable_names);
PrincipalAngleFractionTable = cellRowsToTableLocal( ...
    fraction_rows, fraction_variable_names);
SimilarityTable = cellRowsToTableLocal( ...
    similarity_rows, similarity_variable_names);

SubspaceOverlapSigTest.tables = struct();
SubspaceOverlapSigTest.tables.principal_angles = PrincipalAngleTable;
SubspaceOverlapSigTest.tables.principal_angle_fractions = ...
    PrincipalAngleFractionTable;
SubspaceOverlapSigTest.tables.similarity = SimilarityTable;

output_base = fullfile(result_dir, sprintf( ...
    'subspace_overlap_sigtest_%s_%s', ...
    control_method, latent_selection_tag));

if save_results
    output_mat = [output_base, '.mat'];
    save(output_mat, 'SubspaceOverlapSigTest', '-v7.3');
    fprintf('\nSaved %s\n', output_mat);
end

if save_tables
    principal_csv = [output_base, '_principal_angles.csv'];
    fraction_csv = [output_base, '_principal_angle_fractions.csv'];
    similarity_csv = [output_base, '_similarity.csv'];

    writetable(PrincipalAngleTable, principal_csv);
    writetable(PrincipalAngleFractionTable, fraction_csv);
    writetable(SimilarityTable, similarity_csv);

    fprintf('Saved %s\n', principal_csv);
    fprintf('Saved %s\n', fraction_csv);
    fprintf('Saved %s\n', similarity_csv);
end

fprintf('\nsubspace_overlap_sigtest completed.\n');

%% ======================================================================
% Local functions
%% ======================================================================

function method = normalizeControlMethodLocal(method_in)

method = lower(strtrim(char(string(method_in))));
method_key = regexprep(method, '[^a-z0-9]+', '');

switch method_key
    case {'uniform', 'haar', 'uniformhaar'}
        method = 'uniform';
    case {'covarianceweighted', 'covariance', 'elsayed'}
        method = 'covariance_weighted';
    otherwise
        error(['Unknown control_method: %s. Use ''uniform'' or ', ...
               '''covariance_weighted''.'], char(string(method_in)));
end

end

function [covariance_info, sampling_factor] = loadGroupCovariancesLocal( ...
        dat_file, stim_tag, data_content, group_p, num_groups)

if ~isfile(dat_file) && isfile([dat_file, '.mat'])
    dat_file = [dat_file, '.mat'];
end
if ~isfile(dat_file)
    error('Covariance source file not found: %s', dat_file);
end

fprintf('Loading covariance source: %s\n', dat_file);
Sdata = load(dat_file, 'model_data_allruns');
if ~isfield(Sdata, 'model_data_allruns')
    error('%s does not contain model_data_allruns.', dat_file);
end

model_data_allruns = Sdata.model_data_allruns;
all_tags = getAllRunTagsLocal(model_data_allruns);
run_match = find(strcmp(all_tags, stim_tag));
if isempty(run_match)
    error('stim_tag not found in model_data_allruns: %s', stim_tag);
end
if numel(run_match) > 1
    error('stim_tag occurs more than once in model_data_allruns: %s', stim_tag);
end

this_run = getRunEntryLocal(model_data_allruns, run_match);
if ~isfield(this_run, data_content)
    error('Selected model_data_allruns entry is missing field %s.', data_content);
end

seq_true = this_run.(data_content);
Y = concatenateSequenceDataLocal(seq_true);
if isempty(Y) || any(~isfinite(Y(:)))
    error(['The pooled %s data are empty or contain nonfinite values. ', ...
           'Covariance control requires finite model input data.'], data_content);
end

y_dims = getRunGroupDimsLocal(this_run, data_content);
if numel(y_dims) ~= num_groups
    error('model_data_allruns has %d groups; group_names has %d.', ...
        numel(y_dims), num_groups);
end
if size(Y, 1) ~= sum(y_dims)
    error('Pooled data have %d rows; sum(yDims) is %d.', ...
        size(Y, 1), sum(y_dims));
end
if any(y_dims ~= group_p)
    error(['Neural dimensions in model_data_allruns [%s] do not match ', ...
           'the saved subspace bases [%s].'], ...
        num2str(y_dims), num2str(group_p));
end

obs_start = cumsum([1, y_dims(1:end-1)]);
obs_end = cumsum(y_dims);

covariance_info = cell(1, num_groups);
sampling_factor = cell(1, num_groups);

for g = 1:num_groups
    Yg = double(Y(obs_start(g):obs_end(g), :));
    Sigma = cov(Yg.', 0);
    Sigma = real((Sigma + Sigma.') / 2);

    [V, D] = eig(Sigma);
    eigenvalues = real(diag(D));
    scale = max(1, max(abs(eigenvalues)));
    eigenvalues(eigenvalues < 0 & eigenvalues > -1e-10 * scale) = 0;
    if any(eigenvalues < 0)
        error('Group %d covariance has materially negative eigenvalues.', g);
    end

    sampling_factor{g} = V * diag(sqrt(eigenvalues));
    covariance_info{g} = struct( ...
        'used', true, ...
        'matrix', Sigma, ...
        'num_samples', size(Yg, 2), ...
        'source_file', dat_file, ...
        'stim_tag', stim_tag, ...
        'data_content', data_content, ...
        'sample_definition', 'all trials x all time bins');

    fprintf('  Group %d covariance: p=%d, samples=%d\n', ...
        g, size(Yg, 1), size(Yg, 2));
end

end

function [QA, QB] = drawMatchedSubspacesLocal( ...
        p, k, m, control_method, sampling_factor, max_draw_attempts)

for attempt = 1:max_draw_attempts
    raw_joint = randn(p, k + m);
    if strcmp(control_method, 'covariance_weighted')
        raw_joint = sampling_factor * raw_joint;
    end

    raw_A = raw_joint(:, 1:k);
    raw_B = raw_joint(:, k + (1:m));

    if rank(raw_A) == k && rank(raw_B) == m && ...
            rank(raw_joint) == k + m
        QA = fullColumnBasisLocal(raw_A, k);
        QB = fullColumnBasisLocal(raw_B, m);
        return;
    end
end

error(['Failed to draw individually and jointly full-rank random vectors ', ...
       'after %d attempts (p=%d, k=%d, m=%d).'], ...
    max_draw_attempts, p, k, m);

end

function Q = fullColumnBasisLocal(X, expected_rank)

[U, ~, ~] = svd(X, 'econ');
if size(U, 2) < expected_rank
    error('SVD returned fewer than %d basis vectors.', expected_rank);
end
Q = U(:, 1:expected_rank);

end

function [angle_deg, sim_A_captures_B, sim_B_captures_A] = ...
        calculateAllOverlapMetricsLocal(QA, QB)

principal_cosine = svd(QA.' * QB, 'econ');
principal_cosine = min(max(real(principal_cosine), 0), 1);
angle_deg = reshape(acos(principal_cosine) * 180 / pi, 1, []);

sim_A_captures_B = directionalSimilarityLocal(QA, QB);
sim_B_captures_A = directionalSimilarityLocal(QB, QA);

end

function value = directionalSimilarityLocal(Qcapture, Qtarget)

residual = Qtarget - Qcapture * (Qcapture.' * Qtarget);
value = 1 - norm(residual, 'fro') / norm(Qtarget, 'fro');
value = min(max(real(value), 0), 1);

end

function result = principalAngleTestLocal(observed, null_angles, alpha)

B = size(null_angles, 1);
r = size(null_angles, 2);

result = struct();
result.observed_angle_deg = observed;
result.null_median_deg = median(null_angles, 1);
result.null_ci95_deg = nan(r, 2);
result.p_small = nan(1, r);
result.p_large = nan(1, r);
result.p_two = nan(1, r);
result.direction = cell(1, r);
result.significant_two_sided = false(1, r);
result.significant_aligned = false(1, r);
result.significant_orthogonal = false(1, r);

for j = 1:r
    null_j = null_angles(:, j);
    result.null_ci95_deg(j, :) = empiricalIntervalLocal(null_j, [0.025, 0.975]);
    result.p_small(j) = (1 + sum(null_j <= observed(j))) / (B + 1);
    result.p_large(j) = (1 + sum(null_j >= observed(j))) / (B + 1);
    result.p_two(j) = min(1, 2 * min(result.p_small(j), result.p_large(j)));

    if observed(j) < result.null_median_deg(j)
        result.direction{j} = 'aligned';
    elseif observed(j) > result.null_median_deg(j)
        result.direction{j} = 'orthogonal';
    else
        result.direction{j} = 'at_null_median';
    end

    result.significant_two_sided(j) = result.p_two(j) < alpha;
    result.significant_aligned(j) = ...
        result.significant_two_sided(j) && strcmp(result.direction{j}, 'aligned');
    result.significant_orthogonal(j) = ...
        result.significant_two_sided(j) && strcmp(result.direction{j}, 'orthogonal');
end

result.num_aligned = sum(result.significant_aligned);
result.num_orthogonal = sum(result.significant_orthogonal);
result.fraction_aligned = result.num_aligned / r;
result.fraction_orthogonal = result.num_orthogonal / r;

[result.null_fraction_aligned, result.null_fraction_orthogonal] = ...
    leaveOneOutNullAngleFractionsLocal(null_angles, alpha);

end

function [fraction_aligned, fraction_orthogonal] = ...
        leaveOneOutNullAngleFractionsLocal(null_angles, alpha)

% Each null row is tested against the other B-1 rows. For a value x_b,
% the +1 Monte Carlo numerator after excluding itself equals the count in
% the full B-row sample, so p_small=count(values<=x_b)/B and likewise for
% p_large. This creates one whole-spectrum fraction per matched null draw.

[B, r] = size(null_angles);
aligned = false(B, r);
orthogonal = false(B, r);

for j = 1:r
    values = null_angles(:, j);
    [sorted_values, order] = sort(values, 'ascend');

    count_le_sorted = nan(B, 1);
    count_ge_sorted = nan(B, 1);

    first_idx = 1;
    while first_idx <= B
        last_idx = first_idx;
        while last_idx < B && sorted_values(last_idx + 1) == sorted_values(first_idx)
            last_idx = last_idx + 1;
        end
        count_le_sorted(first_idx:last_idx) = last_idx;
        count_ge_sorted(first_idx:last_idx) = B - first_idx + 1;
        first_idx = last_idx + 1;
    end

    count_le = nan(B, 1);
    count_ge = nan(B, 1);
    count_le(order) = count_le_sorted;
    count_ge(order) = count_ge_sorted;

    p_small = count_le / B;
    p_large = count_ge / B;
    p_two = min(1, 2 * min(p_small, p_large));

    aligned(:, j) = p_two < alpha & p_small < p_large;
    orthogonal(:, j) = p_two < alpha & p_large < p_small;
end

fraction_aligned = mean(aligned, 2);
fraction_orthogonal = mean(orthogonal, 2);

end

function result = scalarControlTestLocal(observed, null_values, alpha)

null_values = null_values(:);
B = numel(null_values);

result = struct();
result.observed = observed;
result.null_median = median(null_values);
result.null_mean = mean(null_values);
result.null_ci95 = empiricalIntervalLocal(null_values, [0.025, 0.975]);
result.effect_from_null_median = observed - result.null_median;
result.p_high = (1 + sum(null_values >= observed)) / (B + 1);
result.p_low = (1 + sum(null_values <= observed)) / (B + 1);
result.p_two = min(1, 2 * min(result.p_high, result.p_low));

if result.effect_from_null_median > 0
    result.direction = 'higher_than_control';
    result.p_directional = result.p_high;
elseif result.effect_from_null_median < 0
    result.direction = 'lower_than_control';
    result.p_directional = result.p_low;
else
    result.direction = 'at_control_median';
    result.p_directional = min(result.p_high, result.p_low);
end

% Per the requested reporting rule, similarity significance and stars use
% the one-sided p-value in the observed direction.
result.significant_directional = result.p_directional < alpha;
result.star = starFromPValueLocal(result.p_directional);

end

function interval = empiricalIntervalLocal(values, probabilities)

values = sort(values(isfinite(values)));
if isempty(values)
    interval = nan(size(probabilities));
    return;
end

n = numel(values);
interval = nan(size(probabilities));

for q = 1:numel(probabilities)
    position = 1 + (n - 1) * probabilities(q);
    lo = floor(position);
    hi = ceil(position);
    if lo == hi
        interval(q) = values(lo);
    else
        weight = position - lo;
        interval(q) = values(lo) * (1 - weight) + values(hi) * weight;
    end
end

end

function star = starFromPValueLocal(p)

if ~isfinite(p)
    star = '';
elseif p < 0.001
    star = '***';
elseif p < 0.01
    star = '**';
elseif p < 0.05
    star = '*';
else
    star = '';
end

end

function pair_out = initializePairOutputLocal(observed_pair, pair_name)

pair_out = struct();
pair_out.name = pair_name;
pair_out.labelA = getTextFieldLocal(observed_pair, 'labelA', 'A');
pair_out.labelB = getTextFieldLocal(observed_pair, 'labelB', 'B');
pair_out.status = getTextFieldLocal(observed_pair, 'status', 'ok');
pair_out.warning = getTextFieldLocal(observed_pair, 'warning', '');
pair_out.p = NaN;
pair_out.k = 0;
pair_out.m = 0;
pair_out.r = 0;
pair_out.capture_field_A_captures_B = '';
pair_out.capture_field_B_captures_A = '';
pair_out.observed = struct();
pair_out.null = struct();
pair_out.principal = struct();
pair_out.similarity = struct();

if ~isfield(observed_pair, 'basisA') || ~isfield(observed_pair, 'basisB') || ...
        isempty(observed_pair.basisA) || isempty(observed_pair.basisB)
    pair_out.status = 'skipped_empty_subspace';
    if isempty(pair_out.warning)
        pair_out.warning = 'Observed pair has an empty subspace.';
    end
end

end

function value = getObservedSimilarityLocal(observed_pair, field_name)

if ~isfield(observed_pair, 'similarity') || ...
        ~isfield(observed_pair.similarity, field_name)
    error('Observed pair is missing similarity.%s.', field_name);
end
value = observed_pair.similarity.(field_name);
if ~(isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value))
    error('Observed similarity.%s must be one finite numeric scalar.', field_name);
end
value = double(value);

end

function row = makeSimilarityTableRowLocal( ...
        group_index, group_name, pair_name, direction_name, k, m, test)

row = { ...
    group_index, group_name, pair_name, direction_name, k, m, ...
    test.observed, test.null_median, test.effect_from_null_median, ...
    test.p_high, test.p_low, test.p_two, test.p_directional, ...
    test.direction, test.significant_directional, test.star};

end

function T = cellRowsToTableLocal(rows, variable_names)

if isempty(rows)
    T = cell2table(cell(0, numel(variable_names)), ...
        'VariableNames', variable_names);
else
    T = cell2table(rows, 'VariableNames', variable_names);
end

end

function p = inferGroupObservationDimLocal(SubspaceSim, group_index)

group_result = SubspaceSim.group(group_index);
p_candidates = [];

for pair_idx = 1:numel(group_result.pair)
    pair_result = group_result.pair{pair_idx};
    basis_fields = {'basisA', 'basisB'};
    for f = 1:numel(basis_fields)
        if isfield(pair_result, basis_fields{f}) && ...
                ~isempty(pair_result.(basis_fields{f}))
            p_candidates(end + 1) = size(pair_result.(basis_fields{f}), 1); %#ok<AGROW>
        end
    end
end

if isempty(p_candidates)
    error('Could not infer neural dimension p for group %d.', group_index);
end
if any(p_candidates ~= p_candidates(1))
    error('Saved bases disagree on neural dimension p for group %d.', group_index);
end
p = p_candidates(1);

end

function validateSubspaceSimLocal(SubspaceSim, num_groups)

if ~isstruct(SubspaceSim) || ~isfield(SubspaceSim, 'group')
    error('SubspaceSim must be a struct containing group.');
end
if numel(SubspaceSim.group) ~= num_groups
    error(['SubspaceSim contains %d groups, but group_names contains %d. ', ...
           'Their order and count must match.'], ...
        numel(SubspaceSim.group), num_groups);
end

required_pairs = {'across_vs_within', 'feedforward_vs_feedback'};
for g = 1:num_groups
    for pair_idx = 1:numel(required_pairs)
        getPairByNameLocal(SubspaceSim, g, required_pairs{pair_idx});
    end
end

end

function pair_result = getPairByNameLocal( ...
        SubspaceSim, group_index, requested_name)

group_result = SubspaceSim.group(group_index);
if ~isfield(group_result, 'pair') || ~iscell(group_result.pair)
    error('SubspaceSim.group(%d).pair must be a cell array.', group_index);
end

pair_index = [];
if isfield(group_result, 'pairNames')
    saved_names = cellstr(string(group_result.pairNames));
    pair_index = find(strcmp(saved_names, requested_name), 1);
end

if isempty(pair_index)
    for p = 1:numel(group_result.pair)
        this_pair = group_result.pair{p};
        if isstruct(this_pair) && isfield(this_pair, 'name') && ...
                strcmp(char(string(this_pair.name)), requested_name)
            pair_index = p;
            break;
        end
    end
end

if isempty(pair_index)
    error('Group %d is missing pair %s.', group_index, requested_name);
end
pair_result = group_result.pair{pair_index};
if ~isstruct(pair_result)
    error('SubspaceSim.group(%d).pair{%d} must be a struct.', ...
        group_index, pair_index);
end

end

function field_name = makeCaptureFieldNameLocal(label_A, label_B)

field_name = sprintf('%s_captures_%s', ...
    normalizeCaptureLabelLocal(label_A), ...
    normalizeCaptureLabelLocal(label_B));

end

function label = normalizeCaptureLabelLocal(label_in)

label = lower(char(string(label_in)));
label = regexprep(label, '[^a-z0-9]+', '_');
label = regexprep(label, '^_+|_+$', '');
if isempty(label)
    label = 'subspace';
end

end

function value = getTextFieldLocal(S, field_name, default_value)

if isfield(S, field_name)
    value = char(string(S.(field_name)));
else
    value = default_value;
end

end

function run_entry = getRunEntryLocal(model_data_allruns, index)

if iscell(model_data_allruns)
    run_entry = model_data_allruns{index};
else
    run_entry = model_data_allruns(index);
end

end

function tags = getAllRunTagsLocal(model_data_allruns)

tags = cell(numel(model_data_allruns), 1);
for j = 1:numel(model_data_allruns)
    run_entry = getRunEntryLocal(model_data_allruns, j);
    if ~isfield(run_entry, 'stim_tag')
        error('stim_tag is missing from model_data_allruns entry %d.', j);
    end
    tags{j} = char(string(run_entry.stim_tag));
end

end

function y_dims = getRunGroupDimsLocal(this_run, data_content)

if isfield(this_run, 'nan_trial_strategy') && ...
        isequal(this_run.nan_trial_strategy, 6)
    field_name = sprintf('%s_groupd', data_content);
    if ~isfield(this_run, field_name)
        error('nan_trial_strategy=6, but %s is missing.', field_name);
    end
    y_dims = this_run.(field_name);
else
    if ~isfield(this_run, 'groupd')
        error('Selected model_data_allruns entry is missing groupd.');
    end
    y_dims = this_run.groupd;
end

y_dims = reshape(double(y_dims), 1, []);
if any(~isfinite(y_dims)) || any(y_dims <= 0) || any(y_dims ~= round(y_dims))
    error('Group dimensions must be finite positive integers.');
end

end

function Y = concatenateSequenceDataLocal(seq_true)

if isempty(seq_true)
    Y = [];
    return;
end

if iscell(seq_true)
    sequence_cells = seq_true(:);
else
    sequence_cells = num2cell(seq_true(:));
end

Y_cells = cell(numel(sequence_cells), 1);
num_rows = [];

for j = 1:numel(sequence_cells)
    this_sequence = sequence_cells{j};
    if ~isstruct(this_sequence) || ~isfield(this_sequence, 'y')
        error('Every sequence must be a struct containing y.');
    end
    this_Y = this_sequence.y;
    if ~(isnumeric(this_Y) && isreal(this_Y) && ismatrix(this_Y))
        error('Every sequence y must be a real numeric matrix.');
    end
    if isempty(num_rows)
        num_rows = size(this_Y, 1);
    elseif size(this_Y, 1) ~= num_rows
        error('Sequence y row counts are inconsistent.');
    end
    Y_cells{j} = double(this_Y);
end

Y = cat(2, Y_cells{:});

end

function names = normalizeGroupNamesLocal(names_in)

if isstring(names_in)
    names = cellstr(names_in(:));
elseif ischar(names_in)
    names = cellstr(names_in);
elseif iscell(names_in)
    names = cell(numel(names_in), 1);
    for k = 1:numel(names_in)
        names{k} = char(string(names_in{k}));
    end
else
    error('group_names must be a char array, string array, or cell array.');
end

names = reshape(names, 1, []);
if isempty(names) || any(cellfun(@isempty, names))
    error('group_names must contain at least one nonempty name.');
end

end

function [display_names, file_tags] = buildGroupLabelsLocal(group_names)

num_groups = numel(group_names);
display_names = cell(1, num_groups);
file_tags = cell(1, num_groups);

for g = 1:num_groups
    display_names{g} = sprintf('Group %d: %s', g, group_names{g});
    file_tags{g} = sprintf('G%02d_%s', g, sanitizeTagLocal(group_names{g}));
end

end

function tag = makeLatentSelectionTagLocal( ...
        use_dsl_filter, dsl_field, ...
        use_svexp_filter, svexp_field, shared_varexp_threshold)

parts = {};
if use_dsl_filter
    parts{end + 1} = sprintf( ...
        'DSL_%s_filtered', sanitizeTagLocal(dsl_field)); %#ok<AGROW>
end
if use_svexp_filter
    parts{end + 1} = sprintf( ...
        'SVExp_%s_%s_filtered', ...
        thresholdTagLocal(shared_varexp_threshold), ...
        sanitizeTagLocal(svexp_field)); %#ok<AGROW>
end
if isempty(parts)
    tag = 'all_latents';
else
    tag = strjoin(parts, '_and_');
end

end

function tag = sanitizeTagLocal(value)

tag = char(string(value));
tag = regexprep(tag, '[^A-Za-z0-9]+', '_');
tag = regexprep(tag, '^_+|_+$', '');
if isempty(tag)
    tag = 'field';
end

end

function tag = thresholdTagLocal(threshold)

tag = sprintf('threshold%.6g', threshold);
tag = strrep(tag, '.', 'p');
tag = strrep(tag, '-', 'm');

end
