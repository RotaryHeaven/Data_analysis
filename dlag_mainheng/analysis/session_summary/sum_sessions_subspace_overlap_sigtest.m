%% sum_sessions_subspace_overlap_sigtest.m
% Across-session matched-control tests for DLAG subspace overlap.
%
% Prerequisite for every session:
%   1) run subspace_similarity_dlag.m;
%   2) run subspace_overlap_sigtest.m with the same settings below.
%
% This summary is intentionally restricted to the pooled all-condition
% model. Each session is one equally weighted experimental observation.
%
% Principal angles:
%   The directional hypothesis is that the observed subspaces are more
%   aligned than control. Each session therefore contributes the fraction
%   of its ordered principal angles that are significantly smaller than
%   their corresponding null angles using a one-sided lower-tail test.
%   This accommodates different r=min(k,m) across sessions. In every
%   group-null iteration, one complete null spectrum row is sampled from
%   each session, its aligned fraction is used, and session fractions are
%   averaged. Thus angles from the same spectrum are never sampled
%   independently.
%
% Directional similarity:
%   The four directions are tested and plotted separately:
%     Across captures Within
%     Within captures Across
%     FF captures FB
%     FB captures FF
%   No average of opposite directions is calculated. The two directions
%   nevertheless use the same per-session null-row indices, so a paired
%   directional average can be derived later from the saved results.
%
% Figures use:
%   colored circles       = observed session values;
%   colored horizontal bar = arithmetic mean across sessions;
%   gray half violin      = empirical group-null distribution;
%   gray vertical line    = group-null 95% range;
%   gray horizontal line  = group-null median.
%
% Only significant results receive a star and p label. Both the aligned-
% fraction plot and all directional-similarity comparisons use the
% prespecified one-sided alternative "higher/more aligned than control"
% and therefore use p_high.

clc;
clear;

%% ========================== USER SETTINGS ==============================

root_dir = 'I:\np_data';

data_content = 'raw_count';
runIdx = 1;

% Display/file labels only. Order must match the saved neural-group order.
group_names = {'V1', 'MT'};

% Must match subspace_similarity_dlag.m and subspace_overlap_sigtest.m.
use_dsl_filter = false;
dsl_field = 'logical';             % usually 'rawlogical' or 'logical'

use_svexp_filter = false;
svexp_field = 'logical';           % 'rawlogical' or 'logical'
shared_varexp_threshold = 0.95;

% Must match the single-session result to be summarized:
%   'uniform'
%   'covariance_weighted'
control_method = 'uniform';

alpha = 0.05;
num_group_null = 10000;
group_random_seed = 2;

% Saving switches.
save_mat = true;
save_tables = true;
save_fig = false;
save_svg = false;
save_png = true;
png_dpi = 300;

figure_visible = 'on';
close_after_save = false;

% First two rows preserve the established Group 1 / Group 2 colors.
group_colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980
];

within_group_spacing = 1.0;
between_group_spacing = 2.0;

% Automatically zoom each y-axis to the range occupied by observed and
% group-null values while respecting the valid boundaries [0,1] and
% [0,100]. If false, the two manual limits below are used.
auto_y_limits = true;
similarity_y_limits = [0, 1];
principal_fraction_y_limits = [0, 100];
y_axis_padding_fraction = 0.12;
minimum_similarity_y_span = 0.10;
minimum_principal_y_span = 10;

similarity_figure_size_inches = [10.8, 4.6];
principal_figure_size_inches = [6.8, 4.6];

point_size = 24;
point_jitter_width = 0.10;
observed_mean_half_width = 0.16;
observed_mean_line_width = 2.2;

null_violin_width = 0.25;
null_violin_offset = 0.06;
null_violin_color = [0.72, 0.72, 0.72];
null_violin_alpha = 0.55;
null_summary_color = [0.35, 0.35, 0.35];

% Independent display switches for summaries drawn over each gray null
% distribution. These options do not change any calculation or saved p.
show_null_median = true;
show_null_ci95 = true;

font_name = 'Arial';
font_size = 9;

%% ======================== VALIDATE SETTINGS ============================

if ~isfolder(root_dir)
    error('root_dir does not exist: %s', root_dir);
end

group_names = normalizeGroupNamesLocal(group_names);
[group_display_names, ~, group_mapping_tag] = ...
    buildGroupLabelsLocal(group_names);
num_groups = numel(group_names);
group_colors = ensureColorRowsLocal(group_colors, num_groups);

control_method = normalizeControlMethodLocal(control_method);

validateattributes(runIdx, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'runIdx');
validateattributes(alpha, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1}, mfilename, 'alpha');
validateattributes(num_group_null, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'num_group_null');
validateattributes(group_random_seed, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', 'finite'}, ...
    mfilename, 'group_random_seed');
validateattributes(png_dpi, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'png_dpi');
validateattributes(auto_y_limits, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'auto_y_limits');
validateattributes(show_null_median, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'show_null_median');
validateattributes(show_null_ci95, {'logical', 'numeric'}, ...
    {'scalar'}, mfilename, 'show_null_ci95');
validateattributes(y_axis_padding_fraction, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnegative'}, ...
    mfilename, 'y_axis_padding_fraction');
validateattributes(minimum_similarity_y_span, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<=', 1}, ...
    mfilename, 'minimum_similarity_y_span');
validateattributes(minimum_principal_y_span, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'positive', '<=', 100}, ...
    mfilename, 'minimum_principal_y_span');

auto_y_limits = logical(auto_y_limits);
show_null_median = logical(show_null_median);
show_null_ci95 = logical(show_null_ci95);

validateManualLimitsLocal(similarity_y_limits, [0, 1], ...
    'similarity_y_limits');
validateManualLimitsLocal(principal_fraction_y_limits, [0, 100], ...
    'principal_fraction_y_limits');

latent_selection_tag = makeLatentSelectionTagLocal( ...
    use_dsl_filter, dsl_field, ...
    use_svexp_filter, svexp_field, shared_varexp_threshold);

output_base = fullfile(root_dir, sprintf( ...
    '%s_all_condition_M_run%03d_subspace_overlap_sigtest_%s_%s_%s', ...
    data_content, runIdx, control_method, ...
    latent_selection_tag, group_mapping_tag));

fprintf('Root dir              : %s\n', root_dir);
fprintf('Model mode            : all_condition_model\n');
fprintf('Control method        : %s\n', control_method);
fprintf('Latent selection tag  : %s\n', latent_selection_tag);
fprintf('Group null iterations : %d\n', num_group_null);
fprintf('Group labels          : %s\n', strjoin(group_display_names, ' | '));
fprintf('Output base           : %s\n', output_base);

%% ======================= FIND / READ SESSIONS ==========================

session_dirs = findCatgtSessionDirsLocal(root_dir);
if isempty(session_dirs)
    error('No catgt_* folders found one level below root_dir: %s', root_dir);
end

fprintf('Found %d candidate catgt_* folders.\n', numel(session_dirs));

records = struct('session_name', {}, 'session_dir', {}, ...
    'input_file', {}, 'group', {});
skipped = struct('session_dir', {}, 'reason', {});

for session_idx = 1:numel(session_dirs)
    session_dir = session_dirs{session_idx};
    [~, session_name] = fileparts(session_dir);
    input_file = makeInputFileLocal( ...
        session_dir, data_content, runIdx, ...
        control_method, latent_selection_tag);

    fprintf('[%d/%d] %s\n', session_idx, numel(session_dirs), session_dir);

    if ~isfile(input_file)
        reason = sprintf('Input file not found: %s', input_file);
        warning('%s', reason);
        skipped(end + 1) = struct( ... %#ok<SAGROW>
            'session_dir', session_dir, 'reason', reason);
        continue;
    end

    try
        S = load(input_file, 'SubspaceOverlapSigTest');
        if ~isfield(S, 'SubspaceOverlapSigTest')
            error('%s does not contain SubspaceOverlapSigTest.', input_file);
        end

        one_result = S.SubspaceOverlapSigTest;
        validateOneSessionResultLocal( ...
            one_result, num_groups, control_method, ...
            latent_selection_tag, alpha);

        record = compactSessionRecordLocal( ...
            one_result, session_name, session_dir, input_file, ...
            num_groups, alpha);
        records(end + 1) = record; %#ok<SAGROW>

    catch ME
        reason = sprintf('%s: %s', ME.identifier, ME.message);
        warning('Skipping %s. %s', session_dir, reason);
        skipped(end + 1) = struct( ... %#ok<SAGROW>
            'session_dir', session_dir, 'reason', reason);
    end
end

if isempty(records)
    error('No valid single-session subspace-overlap test files were found.');
end

fprintf('Included %d session(s); skipped %d.\n', ...
    numel(records), numel(skipped));

%% ========================= GROUP-NULL TESTS ============================

rng(group_random_seed, 'twister');

pair_names = {'across_vs_within', 'feedforward_vs_feedback'};
pair_xlabels = {'Across vs Within', 'FF vs FB'};

Summary = struct();
Summary.meta = struct();
Summary.meta.program = mfilename;
Summary.meta.created = datestr(now, 30);
Summary.meta.root_dir = root_dir;
Summary.meta.model_mode = 'all_condition_model';
Summary.meta.data_content = data_content;
Summary.meta.runIdx = runIdx;
Summary.meta.control_method = control_method;
Summary.meta.alpha = alpha;
Summary.meta.num_group_null = num_group_null;
Summary.meta.group_random_seed = group_random_seed;
Summary.meta.latent_selection_tag = latent_selection_tag;
Summary.meta.group_names = group_names;
Summary.meta.group_display_names = group_display_names;
Summary.meta.group_mapping_tag = group_mapping_tag;
Summary.meta.included_session_count = numel(records);
Summary.meta.group_weighting = 'equal weight per session';
Summary.meta.similarity_reporting = ...
    'one-sided p_high for higher similarity than control; opposite directions kept separate';
Summary.meta.principal_reporting = ...
    ['one-sided aligned-angle fractions (p_small<alpha), followed by ', ...
     'an upper-tail group test of their equal-weight session mean'];
Summary.meta.null_sampling = ...
    ['one complete saved null row per session and pair per iteration; ', ...
     'the same row index is used for principal angles and both similarities'];
Summary.meta.null_density_method = ...
    'bounded kernel-density estimate with reflection boundary correction';
Summary.meta.show_null_median = show_null_median;
Summary.meta.show_null_ci95 = show_null_ci95;
Summary.meta.auto_y_limits = auto_y_limits;
Summary.skipped = skipped;
Summary.sessions = records;

similarity_summary_rows = {};
principal_summary_rows = {};
session_value_rows = {};

for g = 1:num_groups
    group_out = struct();
    group_out.name = group_display_names{g};
    group_out.pairNames = pair_names;
    group_out.pair = cell(1, numel(pair_names));

    for pair_idx = 1:numel(pair_names)
        pair_name = pair_names{pair_idx};

        valid_record_indices = findValidPairRecordsLocal(records, g, pair_idx);
        n_sessions = numel(valid_record_indices);

        pair_out = struct();
        pair_out.name = pair_name;
        pair_out.xlabel = pair_xlabels{pair_idx};
        pair_out.valid_record_indices = valid_record_indices;
        pair_out.n_sessions = n_sessions;

        if n_sessions == 0
            pair_out.status = 'no_valid_sessions';
            pair_out.warning = 'No session had a valid result for this pair.';
            group_out.pair{pair_idx} = pair_out;
            continue;
        end

        pair_out.status = 'ok';
        pair_out.warning = '';

        session_names = cell(n_sessions, 1);
        dims_k = nan(n_sessions, 1);
        dims_m = nan(n_sessions, 1);
        dims_r = nan(n_sessions, 1);
        observed_fraction_aligned = nan(n_sessions, 1);
        observed_similarity_ab = nan(n_sessions, 1);
        observed_similarity_ba = nan(n_sessions, 1);

        null_fraction_aligned = cell(n_sessions, 1);
        null_similarity_ab = cell(n_sessions, 1);
        null_similarity_ba = cell(n_sessions, 1);
        null_draw_indices = zeros(num_group_null, n_sessions, 'uint32');

        label_A = '';
        label_B = '';
        field_ab = '';
        field_ba = '';

        for local_idx = 1:n_sessions
            record_idx = valid_record_indices(local_idx);
            one_pair = records(record_idx).group(g).pair{pair_idx};

            session_names{local_idx} = records(record_idx).session_name;
            dims_k(local_idx) = one_pair.k;
            dims_m(local_idx) = one_pair.m;
            dims_r(local_idx) = one_pair.r;
            observed_fraction_aligned(local_idx) = ...
                one_pair.observed_fraction_aligned;
            observed_similarity_ab(local_idx) = one_pair.observed_similarity_ab;
            observed_similarity_ba(local_idx) = one_pair.observed_similarity_ba;

            null_fraction_aligned{local_idx} = one_pair.null_fraction_aligned;
            null_similarity_ab{local_idx} = one_pair.null_similarity_ab;
            null_similarity_ba{local_idx} = one_pair.null_similarity_ba;

            B_i = numel(one_pair.null_similarity_ab);
            null_draw_indices(:, local_idx) = ...
                uint32(randi(B_i, num_group_null, 1));

            if local_idx == 1
                label_A = one_pair.labelA;
                label_B = one_pair.labelB;
                field_ab = one_pair.field_ab;
                field_ba = one_pair.field_ba;
            end
        end

        group_null_fraction_aligned = meanFromPairedDrawsLocal( ...
            null_fraction_aligned, null_draw_indices);
        group_null_similarity_ab = meanFromPairedDrawsLocal( ...
            null_similarity_ab, null_draw_indices);
        group_null_similarity_ba = meanFromPairedDrawsLocal( ...
            null_similarity_ba, null_draw_indices);

        principal_aligned = upperTailGroupTestLocal( ...
            mean(observed_fraction_aligned), ...
            group_null_fraction_aligned, alpha);
        similarity_ab = upperTailGroupTestLocal( ...
            mean(observed_similarity_ab), group_null_similarity_ab, alpha);
        similarity_ba = upperTailGroupTestLocal( ...
            mean(observed_similarity_ba), group_null_similarity_ba, alpha);

        pair_out.labelA = label_A;
        pair_out.labelB = label_B;
        pair_out.field_A_captures_B = field_ab;
        pair_out.field_B_captures_A = field_ba;
        pair_out.session_names = session_names;
        pair_out.dimA_by_session = dims_k;
        pair_out.dimB_by_session = dims_m;
        pair_out.num_angles_by_session = dims_r;
        pair_out.shared_null_draw_indices = null_draw_indices;

        pair_out.principal = struct();
        pair_out.principal.aligned = principal_aligned;
        pair_out.principal.aligned.session_values = observed_fraction_aligned;
        pair_out.principal.aligned.group_null_values = ...
            group_null_fraction_aligned;

        pair_out.similarity = struct();
        pair_out.similarity.(field_ab) = similarity_ab;
        pair_out.similarity.(field_ab).session_values = observed_similarity_ab;
        pair_out.similarity.(field_ab).group_null_values = ...
            group_null_similarity_ab;
        pair_out.similarity.(field_ba) = similarity_ba;
        pair_out.similarity.(field_ba).session_values = observed_similarity_ba;
        pair_out.similarity.(field_ba).group_null_values = ...
            group_null_similarity_ba;

        similarity_summary_rows(end + 1, :) = ... %#ok<SAGROW>
            makeSimilaritySummaryRowLocal( ...
                g, group_display_names{g}, pair_name, ...
                sprintf('%s captures %s', label_A, label_B), ...
                n_sessions, similarity_ab);
        similarity_summary_rows(end + 1, :) = ... %#ok<SAGROW>
            makeSimilaritySummaryRowLocal( ...
                g, group_display_names{g}, pair_name, ...
                sprintf('%s captures %s', label_B, label_A), ...
                n_sessions, similarity_ba);

        principal_summary_rows(end + 1, :) = ... %#ok<SAGROW>
            makePrincipalSummaryRowLocal( ...
                g, group_display_names{g}, pair_name, ...
                pair_xlabels{pair_idx}, n_sessions, principal_aligned);

        for local_idx = 1:n_sessions
            session_value_rows(end + 1, :) = { ... %#ok<SAGROW>
                session_names{local_idx}, g, group_display_names{g}, ...
                pair_name, dims_k(local_idx), dims_m(local_idx), ...
                dims_r(local_idx), observed_fraction_aligned(local_idx), ...
                observed_similarity_ab(local_idx), ...
                observed_similarity_ba(local_idx)};
        end

        group_out.pair{pair_idx} = pair_out;
    end

    Summary.group(g) = group_out;
end

%% ======================== BUILD OUTPUT TABLES ==========================

similarity_variable_names = { ...
    'GroupIndex', 'GroupName', 'PairName', 'DirectionName', ...
    'NumSessions', 'ObservedMean', 'NullMedian', ...
    'NullCI95Low', 'NullCI95High', 'EffectFromNullMedian', ...
    'PHigh', 'SignificantHigh', 'Star'};

principal_variable_names = { ...
    'GroupIndex', 'GroupName', 'PairName', 'Comparison', ...
    'NumSessions', 'ObservedMeanFraction', 'NullMedianFraction', ...
    'NullCI95Low', 'NullCI95High', 'EffectFromNullMedian', ...
    'PHigh', 'SignificantHigh', 'Star'};

session_variable_names = { ...
    'SessionName', 'GroupIndex', 'GroupName', 'PairName', ...
    'DimA', 'DimB', 'NumAngles', ...
    'FractionAligned', ...
    'SimilarityACapturesB', 'SimilarityBCapturesA'};

SimilarityGroupTestTable = cellRowsToTableLocal( ...
    similarity_summary_rows, similarity_variable_names);
PrincipalAngleGroupTestTable = cellRowsToTableLocal( ...
    principal_summary_rows, principal_variable_names);
SessionValueTable = cellRowsToTableLocal( ...
    session_value_rows, session_variable_names);

Summary.tables = struct();
Summary.tables.similarity_group_tests = SimilarityGroupTestTable;
Summary.tables.principal_angle_group_tests = PrincipalAngleGroupTestTable;
Summary.tables.session_values = SessionValueTable;

similarity_plot_values = collectSimilarityPlotValuesLocal(Summary, num_groups);
principal_plot_values = collectPrincipalPlotValuesLocal(Summary, num_groups);

resolved_similarity_y_limits = resolveYLimitsLocal( ...
    similarity_plot_values, [0, 1], auto_y_limits, ...
    similarity_y_limits, y_axis_padding_fraction, ...
    minimum_similarity_y_span);
resolved_principal_y_limits = resolveYLimitsLocal( ...
    100 * principal_plot_values, [0, 100], auto_y_limits, ...
    principal_fraction_y_limits, y_axis_padding_fraction, ...
    minimum_principal_y_span);

Summary.meta.resolved_similarity_y_limits = resolved_similarity_y_limits;
Summary.meta.resolved_principal_y_limits = resolved_principal_y_limits;

fprintf('Similarity y limits   : [%g, %g]\n', resolved_similarity_y_limits);
fprintf('Principal y limits    : [%g, %g]\n', resolved_principal_y_limits);

%% =============================== PLOTS =================================

similarity_figure = plotSimilaritySummaryLocal( ...
    Summary, num_groups, group_display_names, group_colors, ...
    within_group_spacing, between_group_spacing, ...
    resolved_similarity_y_limits, similarity_figure_size_inches, ...
    point_size, point_jitter_width, ...
    observed_mean_half_width, observed_mean_line_width, ...
    null_violin_width, null_violin_offset, ...
    null_violin_color, null_violin_alpha, null_summary_color, ...
    show_null_median, show_null_ci95, ...
    font_name, font_size, figure_visible);

principal_figure = plotPrincipalSummaryLocal( ...
    Summary, num_groups, group_display_names, group_colors, ...
    within_group_spacing, between_group_spacing, ...
    resolved_principal_y_limits, principal_figure_size_inches, ...
    point_size, point_jitter_width, ...
    observed_mean_half_width, observed_mean_line_width, ...
    null_violin_width, null_violin_offset, ...
    null_violin_color, null_violin_alpha, null_summary_color, ...
    show_null_median, show_null_ci95, ...
    font_name, font_size, figure_visible);

%% ================================ SAVE =================================

if save_mat
    mat_file = [output_base, '.mat'];
    save(mat_file, 'Summary', '-v7.3');
    fprintf('Saved %s\n', mat_file);
end

if save_tables
    similarity_csv = [output_base, '_similarity_group_tests.csv'];
    principal_csv = [output_base, '_principal_angle_group_tests.csv'];
    session_csv = [output_base, '_session_values.csv'];

    writetable(SimilarityGroupTestTable, similarity_csv);
    writetable(PrincipalAngleGroupTestTable, principal_csv);
    writetable(SessionValueTable, session_csv);

    fprintf('Saved %s\n', similarity_csv);
    fprintf('Saved %s\n', principal_csv);
    fprintf('Saved %s\n', session_csv);
end

saveFigureFormatsLocal( ...
    similarity_figure, [output_base, '_similarity'], ...
    save_fig, save_svg, save_png, png_dpi);
saveFigureFormatsLocal( ...
    principal_figure, [output_base, '_principal_angles_aligned'], ...
    save_fig, save_svg, save_png, png_dpi);

if close_after_save
    close(similarity_figure);
    close(principal_figure);
end

fprintf('\nsum_sessions_subspace_overlap_sigtest completed.\n');

%% ======================================================================
% Local functions
%% ======================================================================

function record = compactSessionRecordLocal( ...
        result, session_name, session_dir, input_file, num_groups, alpha)

record = struct();
record.session_name = session_name;
record.session_dir = session_dir;
record.input_file = input_file;
record.group = repmat(struct('pair', {{}}), 1, num_groups);

pair_names = {'across_vs_within', 'feedforward_vs_feedback'};

for g = 1:num_groups
    record.group(g).pair = cell(1, numel(pair_names));
    for pair_idx = 1:numel(pair_names)
        pair_result = getPairByNameLocal(result, g, pair_names{pair_idx});
        compact = struct();
        compact.valid = false;
        compact.status = getTextFieldLocal(pair_result, 'status', 'unknown');
        compact.warning = getTextFieldLocal(pair_result, 'warning', '');
        compact.labelA = getTextFieldLocal(pair_result, 'labelA', 'A');
        compact.labelB = getTextFieldLocal(pair_result, 'labelB', 'B');
        compact.k = getNumericScalarFieldLocal(pair_result, 'k', NaN);
        compact.m = getNumericScalarFieldLocal(pair_result, 'm', NaN);
        compact.r = getNumericScalarFieldLocal(pair_result, 'r', NaN);
        compact.field_ab = getTextFieldLocal( ...
            pair_result, 'capture_field_A_captures_B', '');
        compact.field_ba = getTextFieldLocal( ...
            pair_result, 'capture_field_B_captures_A', '');

        if ~strcmp(compact.status, 'ok')
            record.group(g).pair{pair_idx} = compact;
            continue;
        end

        if isempty(compact.field_ab) || isempty(compact.field_ba)
            error('Group %d pair %s is missing capture field names.', ...
                g, pair_names{pair_idx});
        end
        if ~isfield(pair_result, 'observed') || ...
                ~isfield(pair_result.observed, 'similarity') || ...
                ~isfield(pair_result.observed.similarity, compact.field_ab) || ...
                ~isfield(pair_result.observed.similarity, compact.field_ba)
            error('Group %d pair %s is missing observed similarities.', ...
                g, pair_names{pair_idx});
        end
        if ~isfield(pair_result, 'null') || ...
                ~isfield(pair_result.null, 'similarity') || ...
                ~isfield(pair_result.null.similarity, compact.field_ab) || ...
                ~isfield(pair_result.null.similarity, compact.field_ba)
            error('Group %d pair %s is missing null similarities.', ...
                g, pair_names{pair_idx});
        end

        if ~isfield(pair_result.observed, 'principal_angle_deg') || ...
                ~isfield(pair_result.null, 'principal_angle_deg')
            error(['Group %d pair %s is missing the complete observed ', ...
                   'or null principal-angle spectrum.'], ...
                g, pair_names{pair_idx});
        end

        observed_angles = finiteVectorLocal( ...
            pair_result.observed.principal_angle_deg, ...
            'observed principal-angle spectrum').';
        null_angles = pair_result.null.principal_angle_deg;
        if ~(isnumeric(null_angles) && isreal(null_angles) && ...
                ismatrix(null_angles) && ~isempty(null_angles) && ...
                all(isfinite(null_angles(:))))
            error('Null principal-angle spectra must be one finite matrix.');
        end
        null_angles = double(null_angles);
        if size(null_angles, 2) ~= numel(observed_angles)
            error(['Group %d pair %s has %d observed angles but %d ', ...
                   'columns in its null spectra.'], ...
                g, pair_names{pair_idx}, numel(observed_angles), ...
                size(null_angles, 2));
        end

        [compact.observed_fraction_aligned, ...
            compact.null_fraction_aligned] = ...
            oneSidedAlignedFractionsLocal( ...
                observed_angles, null_angles, alpha);

        compact.observed_similarity_ab = finiteScalarLocal( ...
            pair_result.observed.similarity.(compact.field_ab), ...
            'observed similarity A captures B');
        compact.observed_similarity_ba = finiteScalarLocal( ...
            pair_result.observed.similarity.(compact.field_ba), ...
            'observed similarity B captures A');
        compact.null_similarity_ab = finiteVectorLocal( ...
            pair_result.null.similarity.(compact.field_ab), ...
            'null similarity A captures B');
        compact.null_similarity_ba = finiteVectorLocal( ...
            pair_result.null.similarity.(compact.field_ba), ...
            'null similarity B captures A');

        B = numel(compact.null_similarity_ab);
        if numel(compact.null_similarity_ba) ~= B || ...
                numel(compact.null_fraction_aligned) ~= B
            error(['Group %d pair %s null arrays do not share one ', ...
                   'common draw count.'], g, pair_names{pair_idx});
        end

        compact.valid = true;
        record.group(g).pair{pair_idx} = compact;
    end
end

end

function validateOneSessionResultLocal( ...
        result, num_groups, control_method, latent_selection_tag, alpha)

if ~isstruct(result) || ~isfield(result, 'meta') || ...
        ~isfield(result, 'group')
    error('Single-session result must contain meta and group.');
end
if numel(result.group) ~= num_groups
    error('Single-session result has %d groups; expected %d.', ...
        numel(result.group), num_groups);
end

if ~isfield(result.meta, 'control_method') || ...
        ~strcmp(normalizeControlMethodLocal(result.meta.control_method), ...
                control_method)
    error('Saved control_method does not match %s.', control_method);
end
if ~isfield(result.meta, 'latent_selection_tag') || ...
        ~strcmp(char(string(result.meta.latent_selection_tag)), ...
                latent_selection_tag)
    error('Saved latent-selection tag does not match %s.', ...
        latent_selection_tag);
end
if ~isfield(result.meta, 'alpha') || ...
        abs(double(result.meta.alpha) - alpha) > 10 * eps(max(1, alpha))
    error('Saved alpha does not match %.6g.', alpha);
end

pair_names = {'across_vs_within', 'feedforward_vs_feedback'};
for g = 1:num_groups
    for p = 1:numel(pair_names)
        getPairByNameLocal(result, g, pair_names{p});
    end
end

end

function indices = findValidPairRecordsLocal(records, group_idx, pair_idx)

indices = [];
for s = 1:numel(records)
    one_pair = records(s).group(group_idx).pair{pair_idx};
    if isstruct(one_pair) && isfield(one_pair, 'valid') && one_pair.valid
        indices(end + 1) = s; %#ok<AGROW>
    end
end

end

function group_mean = meanFromPairedDrawsLocal(value_cells, draw_indices)

num_iterations = size(draw_indices, 1);
num_sessions = size(draw_indices, 2);

if numel(value_cells) ~= num_sessions
    error('value_cells and draw_indices have different session counts.');
end

group_mean = zeros(num_iterations, 1);
for s = 1:num_sessions
    values = value_cells{s};
    indices = double(draw_indices(:, s));
    if any(indices < 1) || any(indices > numel(values))
        error('Saved group-null draw index is out of range.');
    end
    group_mean = group_mean + values(indices) / num_sessions;
end

end

function [observed_fraction, null_fractions] = ...
        oneSidedAlignedFractionsLocal( ...
            observed_angles, null_angle_spectra, alpha)

observed_angles = reshape(observed_angles, 1, []);
null_angle_spectra = double(null_angle_spectra);
B = size(null_angle_spectra, 1);
r = size(null_angle_spectra, 2);

if numel(observed_angles) ~= r || B < 1 || r < 1
    error('Observed and null principal-angle spectra are incompatible.');
end

% Observed angle j is aligned when it lies in the prespecified small-angle
% tail of the corresponding ordered null angle j.
observed_p_small = (1 + sum( ...
    null_angle_spectra <= observed_angles, 1)) / (B + 1);
observed_fraction = mean(observed_p_small < alpha);

% Treat each complete null spectrum in turn as a pseudo-observation. For
% angle j, its leave-one-out Monte Carlo p is its upper tied rank divided
% by B. Whole rows are retained, preserving dependence among the r angles.
null_is_aligned = false(B, r);
for j = 1:r
    [sorted_values, sort_order] = sort(null_angle_spectra(:, j));
    upper_tied_rank = zeros(B, 1);
    first_idx = 1;
    while first_idx <= B
        last_idx = first_idx;
        while last_idx < B && ...
                sorted_values(last_idx + 1) == sorted_values(first_idx)
            last_idx = last_idx + 1;
        end
        upper_tied_rank(sort_order(first_idx:last_idx)) = last_idx;
        first_idx = last_idx + 1;
    end
    null_p_small = upper_tied_rank / B;
    null_is_aligned(:, j) = null_p_small < alpha;
end

null_fractions = mean(null_is_aligned, 2);

end

function result = upperTailGroupTestLocal(observed_mean, null_values, alpha)

null_values = null_values(:);
B = numel(null_values);

result = basicGroupTestFieldsLocal(observed_mean, null_values);
result.p_high = (1 + sum(null_values >= observed_mean)) / (B + 1);
result.significant_high = result.p_high < alpha;
result.star = starFromPValueLocal(result.p_high);

end

function values = collectSimilarityPlotValuesLocal(Summary, num_groups)

values = [];
for g = 1:num_groups
    for pair_idx = 1:numel(Summary.group(g).pair)
        pair_result = Summary.group(g).pair{pair_idx};
        if ~isfield(pair_result, 'status') || ...
                ~strcmp(pair_result.status, 'ok')
            continue;
        end
        field_names = { ...
            pair_result.field_A_captures_B, ...
            pair_result.field_B_captures_A};
        for direction_idx = 1:numel(field_names)
            one_test = pair_result.similarity.(field_names{direction_idx});
            values = [values; one_test.session_values(:); ... %#ok<AGROW>
                one_test.group_null_values(:)];
        end
    end
end

values = values(isfinite(values));

end

function values = collectPrincipalPlotValuesLocal(Summary, num_groups)

values = [];
for g = 1:num_groups
    for pair_idx = 1:numel(Summary.group(g).pair)
        pair_result = Summary.group(g).pair{pair_idx};
        if ~isfield(pair_result, 'status') || ...
                ~strcmp(pair_result.status, 'ok')
            continue;
        end
        one_test = pair_result.principal.aligned;
        values = [values; one_test.session_values(:); ... %#ok<AGROW>
            one_test.group_null_values(:)];
    end
end

values = values(isfinite(values));

end

function y_limits = resolveYLimitsLocal( ...
        values, hard_limits, use_auto, manual_limits, ...
        padding_fraction, minimum_span)

if ~use_auto
    y_limits = double(manual_limits(:)).';
    return;
end

values = values(isfinite(values));
if isempty(values)
    y_limits = hard_limits;
    return;
end

lower_bound = hard_limits(1);
upper_bound = hard_limits(2);
tolerance = 1e-10 * max(1, upper_bound - lower_bound);
if any(values < lower_bound - tolerance) || ...
        any(values > upper_bound + tolerance)
    error('Values used for automatic y limits exceed valid boundaries.');
end
values = min(max(values, lower_bound), upper_bound);

data_low = min(values);
data_high = max(values);
data_span = data_high - data_low;
reference_span = max(data_span, minimum_span);
padding = padding_fraction * reference_span;
y_limits = [data_low - padding, data_high + padding];

if diff(y_limits) < minimum_span
    center = (data_low + data_high) / 2;
    y_limits = center + [-0.5, 0.5] * minimum_span;
end

% Shift, rather than simply truncate, a requested span at a boundary.
if y_limits(1) < lower_bound
    y_limits = y_limits + (lower_bound - y_limits(1));
end
if y_limits(2) > upper_bound
    y_limits = y_limits - (y_limits(2) - upper_bound);
end
y_limits(1) = max(y_limits(1), lower_bound);
y_limits(2) = min(y_limits(2), upper_bound);

if y_limits(2) <= y_limits(1)
    y_limits = hard_limits;
end

end

function validateManualLimitsLocal(limits, hard_limits, variable_name)

if ~(isnumeric(limits) && isreal(limits) && numel(limits) == 2 && ...
        all(isfinite(limits(:))))
    error('%s must contain two finite real numeric values.', variable_name);
end
limits = double(limits(:)).';
if limits(1) >= limits(2)
    error('%s must be strictly increasing.', variable_name);
end
if limits(1) < hard_limits(1) || limits(2) > hard_limits(2)
    error('%s must lie within [%g, %g].', ...
        variable_name, hard_limits(1), hard_limits(2));
end

end

function result = basicGroupTestFieldsLocal(observed_mean, null_values)

result = struct();
result.observed_mean = observed_mean;
result.null_median = median(null_values);
result.null_mean = mean(null_values);
result.null_ci95 = empiricalIntervalLocal(null_values, [0.025, 0.975]);
result.effect_from_null_median = observed_mean - result.null_median;

end

function row = makeSimilaritySummaryRowLocal( ...
        group_idx, group_name, pair_name, direction_name, ...
        num_sessions, test)

row = { ...
    group_idx, group_name, pair_name, direction_name, num_sessions, ...
    test.observed_mean, test.null_median, ...
    test.null_ci95(1), test.null_ci95(2), ...
    test.effect_from_null_median, test.p_high, ...
    test.significant_high, test.star};

end

function row = makePrincipalSummaryRowLocal( ...
        group_idx, group_name, pair_name, comparison, num_sessions, test)

row = { ...
    group_idx, group_name, pair_name, comparison, num_sessions, ...
    test.observed_mean, test.null_median, ...
    test.null_ci95(1), test.null_ci95(2), ...
    test.effect_from_null_median, test.p_high, ...
    test.significant_high, test.star};

end

function fig = plotSimilaritySummaryLocal( ...
        Summary, num_groups, group_display_names, group_colors, ...
        within_spacing, between_spacing, y_limits, figure_size, ...
        point_size, jitter_width, mean_half_width, mean_line_width, ...
        violin_width, violin_offset, violin_color, violin_alpha, ...
        null_line_color, show_null_median, show_null_ci95, ...
        font_name, font_size, figure_visible)

[x_positions, x_labels, category_group, category_pair, category_direction] = ...
    buildSimilarityCategoriesLocal( ...
        num_groups, within_spacing, between_spacing);

fig = figure('Color', 'w', 'Visible', figure_visible, ...
    'Units', 'inches', 'Position', [1, 1, figure_size]);
ax = axes(fig);
hold(ax, 'on');

for c = 1:numel(x_positions)
    g = category_group(c);
    pair_idx = category_pair(c);
    direction_idx = category_direction(c);
    pair_result = Summary.group(g).pair{pair_idx};

    if ~isfield(pair_result, 'status') || ~strcmp(pair_result.status, 'ok')
        continue;
    end

    if direction_idx == 1
        field_name = pair_result.field_A_captures_B;
    else
        field_name = pair_result.field_B_captures_A;
    end
    test = pair_result.similarity.(field_name);

    drawNullDistributionLocal( ...
        ax, x_positions(c), test.group_null_values, ...
        violin_width, violin_offset, violin_color, ...
        violin_alpha, null_line_color, [0, 1], ...
        show_null_median, show_null_ci95);
    drawObservedSessionsLocal( ...
        ax, x_positions(c), test.session_values, group_colors(g, :), ...
        point_size, jitter_width, mean_half_width, mean_line_width);

    if test.significant_high
        text_y = significanceTextYLocal( ...
            test.session_values, test.group_null_values, y_limits);
        text(ax, x_positions(c), text_y, sprintf('%s p_{high}=%s', ...
            test.star, formatPValueLocal(test.p_high)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', font_name, 'FontSize', font_size - 1, ...
            'Interpreter', 'tex', 'Clipping', 'on');
    end
end

formatSummaryAxesLocal( ...
    ax, x_positions, x_labels, y_limits, 'Subspace similarity', ...
    font_name, font_size);

legend_handles = gobjects(num_groups + 1, 1);
for g = 1:num_groups
    legend_handles(g) = plot(ax, nan, nan, 'o', ...
        'MarkerSize', 5.5, ...
        'MarkerFaceColor', group_colors(g, :), ...
        'MarkerEdgeColor', group_colors(g, :), ...
        'LineStyle', 'none');
end
legend_handles(end) = patch(ax, nan, nan, violin_color, ...
    'EdgeColor', 'none', 'FaceAlpha', violin_alpha);
legend_labels = [group_display_names(:); {'Matched control'}];
legend(ax, legend_handles, legend_labels, ...
    'Location', 'northeastoutside', 'Box', 'off', ...
    'FontName', font_name, 'FontSize', font_size);

hold(ax, 'off');

end

function fig = plotPrincipalSummaryLocal( ...
        Summary, num_groups, group_display_names, group_colors, ...
        within_spacing, between_spacing, y_limits, figure_size, ...
        point_size, jitter_width, mean_half_width, mean_line_width, ...
        violin_width, violin_offset, violin_color, violin_alpha, ...
        null_line_color, show_null_median, show_null_ci95, ...
        font_name, font_size, figure_visible)

[x_positions, x_labels, category_group, category_pair] = ...
    buildPrincipalCategoriesLocal( ...
        num_groups, within_spacing, between_spacing);

fig = figure('Color', 'w', 'Visible', figure_visible, ...
    'Units', 'inches', 'Position', [1, 1, figure_size]);

ax = axes(fig);
hold(ax, 'on');

for c = 1:numel(x_positions)
    g = category_group(c);
    pair_idx = category_pair(c);
    pair_result = Summary.group(g).pair{pair_idx};

    if ~isfield(pair_result, 'status') || ~strcmp(pair_result.status, 'ok')
        continue;
    end

    test = pair_result.principal.aligned;
    observed_percent = 100 * test.session_values;
    null_percent = 100 * test.group_null_values;

    drawNullDistributionLocal( ...
        ax, x_positions(c), null_percent, ...
        violin_width, violin_offset, violin_color, ...
        violin_alpha, null_line_color, [0, 100], ...
        show_null_median, show_null_ci95);
    drawObservedSessionsLocal( ...
        ax, x_positions(c), observed_percent, group_colors(g, :), ...
        point_size, jitter_width, mean_half_width, mean_line_width);

    if test.significant_high
        text_y = significanceTextYLocal( ...
            observed_percent, null_percent, y_limits);
        text(ax, x_positions(c), text_y, sprintf('%s p_{high}=%s', ...
            test.star, formatPValueLocal(test.p_high)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', font_name, 'FontSize', font_size - 1, ...
            'Interpreter', 'tex', 'Clipping', 'on');
    end
end

formatSummaryAxesLocal( ...
    ax, x_positions, x_labels, y_limits, ...
    'Significantly aligned ordered principal angles (%)', ...
    font_name, font_size);
title(ax, 'More aligned than matched control', ...
    'FontName', font_name, 'FontSize', font_size + 1, ...
    'FontWeight', 'normal');

legend_handles = gobjects(num_groups + 1, 1);
for g = 1:num_groups
    legend_handles(g) = plot(ax, nan, nan, 'o', ...
        'MarkerSize', 5.5, ...
        'MarkerFaceColor', group_colors(g, :), ...
        'MarkerEdgeColor', group_colors(g, :), ...
        'LineStyle', 'none');
end
legend_handles(end) = patch(ax, nan, nan, violin_color, ...
    'EdgeColor', 'none', 'FaceAlpha', violin_alpha);
legend_labels = [group_display_names(:); {'Matched control'}];
legend(ax, legend_handles, legend_labels, ...
    'Location', 'northeastoutside', 'Box', 'off', ...
    'FontName', font_name, 'FontSize', font_size);

hold(ax, 'off');

end

function drawNullDistributionLocal( ...
        ax, x, values, width, offset, face_color, face_alpha, line_color, ...
        support_limits, show_median, show_ci95)

values = values(isfinite(values));
if isempty(values)
    return;
end

[density, grid_values] = estimateBoundedDensityLocal( ...
    values, support_limits);
if isempty(density) || max(density) <= 0
    density = ones(size(grid_values));
end
density = density / max(density);

base_x = x + offset;
x_outer = base_x + width * density;

patch(ax, [base_x * ones(size(grid_values)), fliplr(x_outer)], ...
    [grid_values, fliplr(grid_values)], face_color, ...
    'EdgeColor', 'none', 'FaceAlpha', face_alpha, ...
    'HandleVisibility', 'off');

summary_x = base_x + 0.55 * width;

if show_ci95
    ci95 = empiricalIntervalLocal(values, [0.025, 0.975]);
    plot(ax, [summary_x, summary_x], ci95, '-', ...
        'Color', line_color, 'LineWidth', 1.1, ...
        'HandleVisibility', 'off');
end

if show_median
    null_median = median(values);
    plot(ax, [base_x + 0.03 * width, base_x + 0.96 * width], ...
        [null_median, null_median], '-', ...
        'Color', line_color, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

end

function drawObservedSessionsLocal( ...
        ax, x, values, color, point_size, jitter_width, ...
        mean_half_width, mean_line_width)

values = values(isfinite(values));
if isempty(values)
    return;
end

jitter = deterministicJitterLocal(numel(values), jitter_width);
scatter(ax, x - 0.10 + jitter, values, point_size, ...
    'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', color, ...
    'LineWidth', 0.5, 'HandleVisibility', 'off');

observed_mean = mean(values);
plot(ax, [x - mean_half_width - 0.10, x + mean_half_width - 0.10], ...
    [observed_mean, observed_mean], '-', ...
    'Color', color, 'LineWidth', mean_line_width, ...
    'HandleVisibility', 'off');

end

function jitter = deterministicJitterLocal(n, width)

if n <= 1
    jitter = 0;
else
    jitter = linspace(-width / 2, width / 2, n)';
end

end

function [density, grid_values] = estimateBoundedDensityLocal( ...
        values, support_limits)

values = values(:);
lower_bound = support_limits(1);
upper_bound = support_limits(2);
support_span = upper_bound - lower_bound;
tolerance = 1e-10 * max(1, support_span);

if any(values < lower_bound - tolerance) || ...
        any(values > upper_bound + tolerance)
    error('Null values lie outside their declared plotting support.');
end

% Clamp only round-off at an exact mathematical boundary; this never
% changes a substantively out-of-range value because those error above.
values = min(max(values, lower_bound), upper_bound);
grid_values = linspace(lower_bound, upper_bound, 240);

% Silverman's robust bandwidth, bounded away from zero for discrete
% fraction distributions and bounded above to retain visible structure.
n = numel(values);
sample_scale = std(values, 0);
robust_scale = interquartileRangeLocal(values) / 1.349;
positive_scales = [sample_scale, robust_scale];
positive_scales = positive_scales( ...
    isfinite(positive_scales) & positive_scales > 0);
if isempty(positive_scales)
    scale = support_span / 100;
else
    scale = min(positive_scales);
end
bandwidth = 0.9 * scale * n^(-1 / 5);
bandwidth = max(bandwidth, support_span / 500);
bandwidth = min(bandwidth, support_span / 8);

% Reflection correction: for every sample, add its mirror image across
% each finite boundary. The density is evaluated only on the legal
% support, so a distribution of proportions can never appear below 0 or
% above 1 (and a percentage never below 0 or above 100).
density = zeros(size(grid_values));
normalizer = n * bandwidth * sqrt(2 * pi);
chunk_size = 1000;
for first_idx = 1:chunk_size:n
    last_idx = min(n, first_idx + chunk_size - 1);
    one_chunk = values(first_idx:last_idx).';
    z_original = (grid_values(:) - one_chunk) / bandwidth;
    z_lower = (grid_values(:) - (2 * lower_bound - one_chunk)) / ...
        bandwidth;
    z_upper = (grid_values(:) - (2 * upper_bound - one_chunk)) / ...
        bandwidth;
    density = density + sum( ...
        exp(-0.5 * z_original.^2) + ...
        exp(-0.5 * z_lower.^2) + ...
        exp(-0.5 * z_upper.^2), 2).';
end
density = density / normalizer;

end

function formatSummaryAxesLocal( ...
        ax, x_positions, x_labels, y_limits, y_label, font_name, font_size)

set(ax, ...
    'XTick', x_positions, ...
    'XTickLabel', x_labels, ...
    'TickLabelInterpreter', 'none', ...
    'FontName', font_name, ...
    'FontSize', font_size, ...
    'LineWidth', 0.75, ...
    'TickDir', 'out', ...
    'Box', 'off');

xtickangle(ax, 30);
xlim(ax, [min(x_positions) - 0.55, max(x_positions) + 0.55]);
ylim(ax, y_limits);
ylabel(ax, y_label, 'FontName', font_name, 'FontSize', font_size);

end

function y = significanceTextYLocal(observed_values, null_values, y_limits)

all_values = [observed_values(:); null_values(:)];
all_values = all_values(isfinite(all_values));
span = y_limits(2) - y_limits(1);

if isempty(all_values)
    y = y_limits(1) + 0.90 * span;
else
    y = max(all_values) + 0.025 * span;
    y = min(y, y_limits(2) - 0.045 * span);
end

end

function [x_positions, x_labels, category_group, ...
        category_pair, category_direction] = ...
        buildSimilarityCategoriesLocal( ...
            num_groups, within_spacing, between_spacing)

labels_one_group = { ...
    'Across captures Within', ...
    'Within captures Across', ...
    'FF captures FB', ...
    'FB captures FF'};

x_positions = [];
x_labels = {};
category_group = [];
category_pair = [];
category_direction = [];

next_start = 1;
for g = 1:num_groups
    group_x = next_start + (0:3) * within_spacing;
    x_positions = [x_positions, group_x]; %#ok<AGROW>
    x_labels = [x_labels, labels_one_group]; %#ok<AGROW>
    category_group = [category_group, repmat(g, 1, 4)]; %#ok<AGROW>
    category_pair = [category_pair, 1, 1, 2, 2]; %#ok<AGROW>
    category_direction = [category_direction, 1, 2, 1, 2]; %#ok<AGROW>
    next_start = group_x(end) + between_spacing;
end

end

function [x_positions, x_labels, category_group, category_pair] = ...
        buildPrincipalCategoriesLocal( ...
            num_groups, within_spacing, between_spacing)

labels_one_group = {'Across vs Within', 'FF vs FB'};

x_positions = [];
x_labels = {};
category_group = [];
category_pair = [];

next_start = 1;
for g = 1:num_groups
    group_x = next_start + (0:1) * within_spacing;
    x_positions = [x_positions, group_x]; %#ok<AGROW>
    x_labels = [x_labels, labels_one_group]; %#ok<AGROW>
    category_group = [category_group, g, g]; %#ok<AGROW>
    category_pair = [category_pair, 1, 2]; %#ok<AGROW>
    next_start = group_x(end) + between_spacing;
end

end

function saveFigureFormatsLocal( ...
        fig, output_stem, save_fig, save_svg, save_png, png_dpi)

if save_fig
    fig_file = [output_stem, '.fig'];
    savefig(fig, fig_file);
    fprintf('Saved %s\n', fig_file);
end

if save_svg
    svg_file = [output_stem, '.svg'];
    try
        exportgraphics(fig, svg_file, 'ContentType', 'vector');
    catch
        print(fig, svg_file, '-dsvg', '-painters');
    end
    fprintf('Saved %s\n', svg_file);
end

if save_png
    png_file = [output_stem, '.png'];
    try
        exportgraphics(fig, png_file, 'Resolution', png_dpi);
    catch
        print(fig, png_file, '-dpng', sprintf('-r%d', png_dpi));
    end
    fprintf('Saved %s\n', png_file);
end

end

function result = getPairByNameLocal( ...
        saved_result, group_index, requested_name)

group_result = saved_result.group(group_index);
if ~isfield(group_result, 'pair') || ~iscell(group_result.pair)
    error('group(%d).pair must be a cell array.', group_index);
end

pair_index = [];
if isfield(group_result, 'pairNames')
    pair_names = cellstr(string(group_result.pairNames));
    pair_index = find(strcmp(pair_names, requested_name), 1);
end

if isempty(pair_index)
    for p = 1:numel(group_result.pair)
        one_pair = group_result.pair{p};
        if isstruct(one_pair) && isfield(one_pair, 'name') && ...
                strcmp(char(string(one_pair.name)), requested_name)
            pair_index = p;
            break;
        end
    end
end

if isempty(pair_index)
    error('Group %d is missing pair %s.', group_index, requested_name);
end

result = group_result.pair{pair_index};
if ~isstruct(result)
    error('group(%d).pair{%d} must be a struct.', group_index, pair_index);
end

end

function input_file = makeInputFileLocal( ...
        session_dir, data_content, runIdx, control_method, latent_tag)

input_file = fullfile( ...
    session_dir, ...
    ['FA_Dlag_', data_content], ...
    'mat_results', ...
    sprintf('run%03d', runIdx), ...
    sprintf('subspace_overlap_sigtest_%s_%s.mat', ...
        control_method, latent_tag));

end

function session_dirs = findCatgtSessionDirsLocal(root_dir)

listing = dir(fullfile(root_dir, '*', 'catgt_*'));
session_dirs = {};

for k = 1:numel(listing)
    if listing(k).isdir
        session_dirs{end + 1, 1} = fullfile( ...
            listing(k).folder, listing(k).name); %#ok<AGROW>
    end
end

if isempty(session_dirs)
    return;
end

session_dirs = unique(session_dirs, 'stable');
lower_dirs = cellfun(@lower, session_dirs, 'UniformOutput', false);
[~, order] = sort(lower_dirs);
session_dirs = session_dirs(order);

end

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

function value = finiteScalarLocal(value, label)

if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
    error('%s must be one finite numeric scalar.', label);
end
value = double(value);

end

function values = finiteVectorLocal(values, label)

if ~(isnumeric(values) && isreal(values) && isvector(values) && ...
        ~isempty(values) && all(isfinite(values(:))))
    error('%s must be a nonempty finite numeric vector.', label);
end
values = double(values(:));

end

function value = getNumericScalarFieldLocal(S, field_name, default_value)

if isfield(S, field_name)
    value = S.(field_name);
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value))
        error('%s must be one finite numeric scalar.', field_name);
    end
    value = double(value);
else
    value = default_value;
end

end

function value = getTextFieldLocal(S, field_name, default_value)

if isfield(S, field_name)
    value = char(string(S.(field_name)));
else
    value = default_value;
end

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

function value = interquartileRangeLocal(values)

quartiles = empiricalIntervalLocal(values, [0.25, 0.75]);
value = quartiles(2) - quartiles(1);

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

function text_value = formatPValueLocal(p)

if p < 0.001
    text_value = sprintf('%.1e', p);
else
    text_value = sprintf('%.3f', p);
end

end

function T = cellRowsToTableLocal(rows, variable_names)

if isempty(rows)
    T = cell2table(cell(0, numel(variable_names)), ...
        'VariableNames', variable_names);
else
    T = cell2table(rows, 'VariableNames', variable_names);
end

end

function colors = ensureColorRowsLocal(colors, num_groups)

if ~(isnumeric(colors) && size(colors, 2) == 3 && ...
        all(isfinite(colors(:))) && all(colors(:) >= 0) && ...
        all(colors(:) <= 1))
    error('group_colors must be a finite N-by-3 RGB matrix in [0,1].');
end

if size(colors, 1) < num_groups
    extra = lines(num_groups);
    colors = [colors; extra((size(colors, 1) + 1):num_groups, :)];
else
    colors = colors(1:num_groups, :);
end

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

function [display_names, file_tags, mapping_tag] = ...
        buildGroupLabelsLocal(group_names)

num_groups = numel(group_names);
display_names = cell(1, num_groups);
file_tags = cell(1, num_groups);

for g = 1:num_groups
    display_names{g} = sprintf('Group %d: %s', g, group_names{g});
    file_tags{g} = sprintf('G%02d_%s', g, sanitizeTagLocal(group_names{g}));
end
mapping_tag = strjoin(file_tags, '_');

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
