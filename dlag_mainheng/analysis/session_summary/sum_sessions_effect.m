%% sum_sessions_effect.m
% Summarize neuron-wise size or contrast effects across sessions.
% Revision: 2026-09-08a
%
% This script reads the result mats written by plot_size_effect.m or
% plot_contrast_effect.m and produces seven independent figures:
%   1. Original vs Use all
%   2. Use all vs Across
%   3. Use all vs Within
%   4. Use all vs Feedforward
%   5. Use all vs Feedback
%   6. Across/Within/Feedforward/Feedback slope summary
%   7. Across/Within/Feedforward/Feedback delta summary
%
% Scatter figures pool neurons from all sessions in one axis. Group is
% encoded by hue and session by color depth. Panels 2-5 can optionally draw
% one regression line for every session x group. Panel 1 never draws a fit.
%
% Slope is calculated separately for every session x group. Both all-point
% and central regressions are saved; regression_mode selects which slope is
% plotted and which optional fit line is drawn. The displayed fit is clipped
% to the central unbroken square and extended to that square's boundary.
%
% Delta is comparison metric minus Use-all metric. Neurons are separated by
% the sign of the Use-all metric. All neurons from all sessions are pooled
% directly for the displayed mean/median; sessions are not averaged first.
% Positive delta extremes can be collapsed into a labeled top overflow level
% independently for the two base-sign rows. Overflow labels use precision
% chosen from each row's own scale, and no break slashes are drawn. This
% changes display only.

clc;
clear;

%% ========================== USER SETTINGS ==============================

root_dir = 'I:\np_data';

% -------------------------------------------------------------------------
% Effect, data, and model settings
% -------------------------------------------------------------------------

effect_type = 'contrast';
% Options:
%   'size'
%   'contrast'

data_content = 'raw_count';
% Examples:
% raw_count, raw_fr, z_within_trial, z_within_condition,
% z_across_conditions, demean_count_within_trial,
% demean_fr_within_trial, demean_pooledsd_within_condition

model_mode = 'all_condition_model';
% Options:
%   'all_condition_model'
%   'condition_specific_models'

runIdx = 1;

effect_metric = 'auto';
% 'auto' uses:
%   size     -> S_norm_diff
%   contrast -> CMI
%
% Explicit size options:
%   'classic_SI', 'delta_SL', 'S_norm_diff'
% Explicit contrast options:
%   'delta_HL', 'CMI'

% Used only when effect_type = 'size'. It must match plot_size_effect.m.
%   0 = all contrasts
%   1 = low contrast only
%   2 = high contrast only
pick_contrast = 0;

% These must match plot_size_effect.m / plot_contrast_effect.m.
analysis_window = 0.4;  % seconds
sliding_step = 0.2;     % seconds
time_window_index = 1;

% A session is skipped if any one of the six required result files is
% absent or incompatible. Set false to stop immediately on the first error.
skip_invalid_sessions = true;

% Fixed analysis fields requested for this summary.
field_original = 'y';
field_all = 'yRecon_use_all';
field_across = 'yRecon_use_across';
field_within = 'yRecon_use_within';
field_feedforward = 'yRecon_use_feedforward';
field_feedback = 'yRecon_use_feedback';

% -------------------------------------------------------------------------
% Session and group appearance
% -------------------------------------------------------------------------

% Manual area labels in model-group order. These labels are used to select
% the matching area-tagged upstream result files and to label/save the
% summary. They do not determine neuron membership and are not compared with
% any group-name field stored inside the input MAT files.
group_names = {'V1', 'MT'};

% Leave empty to label valid sessions S1, S2, ... in alphabetical folder
% order. Or provide one label for each valid session after skipping.
session_display_names = {};

% The first rows are assigned to the first groups. If nGroups is larger than
% the number of supplied rows, additional distinct colors are generated.
group_colors_base = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980];

% Session 1 is the lightest and the last session is the darkest. The same
% mapping is used in all seven figures.
session_shade_range = [0.48, 1.00];

% -------------------------------------------------------------------------
% Scatter and regression settings
% -------------------------------------------------------------------------

show_session_fit_lines = true;
% false: panels 2-5 show points only.
% true : panels 2-5 show one line per session x group.
% Panel 1 never shows regression lines.

regression_mode = 'all';
% 'all'     : use every finite neuron pair.
% 'central' : use the non-tail points defined by the session-specific
%             broken-axis threshold. If no break is active, central = all.
% This same choice controls panel 6 and the optional lines in panels 2-5.

use_broken_axis_display = false;
break_start_prctile = 93.0;
broken_axis_trigger_ratio = 1.0;
tail_display_frac = 0.08;
break_gap_frac = 0.03;
min_tail_display_slope = 0.001;
force_symmetric_broken_axis = false;
draw_broken_axis_marks = true;

scatter_point_size = 28;
scatter_point_alpha = 0.64;
fit_line_width = 1.35;
reference_line_width = 0.9;
zero_line_width = 0.8;
scatter_figure_size_inches = [6.4, 4.6];
scatter_show_legend = true;

% -------------------------------------------------------------------------
% Slope-summary settings
% -------------------------------------------------------------------------

slope_plot_style = 'points_median';
% Options:
%   'bar_points'    : mean bar + session points + median square
%   'violin'        : violin + session points + optional mean/median
%   'points_median' : session points + optional mean/median, no outline

slope_show_mean = struct( ...
    'bar_points', true, ...
    'violin', false, ...
    'points_median', false);

slope_show_median = struct( ...
    'bar_points', false, ...
    'violin', true, ...
    'points_median', true);

slope_y_limits = [];          % [] = automatic
slope_figure_size_inches = [7.4, 4.3];
slope_bar_width = 0.30;
slope_bar_alpha = 0.88;
slope_point_size = 32;
slope_point_alpha = 0.82;
slope_point_jitter_width = 0.070;
slope_violin_width = 0.145;
slope_violin_alpha = 0.20;
slope_median_line_half_width = 0.12;
slope_mean_marker_size = 7;
slope_group_offset = 0.19;

% -------------------------------------------------------------------------
% Delta-summary settings
% -------------------------------------------------------------------------

delta_plot_style = 'points_median';
% Options:
%   'bar_points'    : mean bar + raw points + median square
%   'violin'        : violin + raw points + median line
%   'points_median' : raw points + median line, no violin outline

% Defaults match the subspace summary programs: the bar version displays
% its mean, while violin and points-only versions do not display a mean.
delta_show_mean = struct( ...
    'bar_points', true, ...
    'violin', false, ...
    'points_median', true);

delta_show_median = struct( ...
    'bar_points', false, ...
    'violin', true, ...
    'points_median', false);

% [] calculates the two sign rows independently. A 1 x 2 vector applies
% the same limits to both rows; a 2 x 2 matrix sets one row per sign group.
% Manual limits are used for rows where overflow is inactive.
delta_y_limits = [];

% Collapse the positive extreme tail into one top overflow level. The two
% sign rows are processed independently. Use NaN for automatic thresholds
% or provide [threshold_base_ge_0, threshold_base_lt_0].
use_delta_overflow_display = true;
delta_overflow_prctile = 93.0;
delta_overflow_trigger_ratio = 1.0;
delta_overflow_thresholds = [NaN, NaN];
delta_overflow_gap_frac = 0.080;
delta_overflow_top_pad_frac = 0.060;

delta_figure_size_inches = [8.0, 6.4];
delta_bar_width = 0.30;
delta_bar_alpha = 0.88;
delta_point_size = 18;
delta_point_alpha = 0.58;
delta_point_jitter_width = 0.105;
delta_violin_width = 0.145;
delta_violin_alpha = 0.20;
delta_median_line_half_width = 0.12;
delta_mean_marker_size = 7;
delta_group_offset = 0.19;

% -------------------------------------------------------------------------
% General figure and save settings
% -------------------------------------------------------------------------

figure_visible = 'on';
close_after_save = true;

% Each figure format can be enabled or disabled independently.
save_fig = false;
save_svg = false;
save_png = true;
png_dpi = 300;

save_mat = true;

font_name = 'Arial';
font_size = 9;
axis_line_width = 0.9;

% Within Across/Within and within Feedforward/Feedback, categories are
% close. A larger gap separates the two conceptual levels.
comparison_x_positions = [1.0, 2.0, 3.5, 4.5];

%% ======================== VALIDATE SETTINGS ============================

if ~isfolder(root_dir)
    error('root_dir does not exist: %s', root_dir);
end

effect_type = normalizeEffectTypeLocal(effect_type);
model_mode = normalizeModelModeLocal(model_mode);
effect_metric = resolveEffectMetricLocal(effect_type, effect_metric);
data_content = char(string(data_content));
slope_plot_style = lower(strtrim(char(slope_plot_style)));
delta_plot_style = lower(strtrim(char(delta_plot_style)));
regression_mode = lower(strtrim(char(regression_mode)));

if isempty(strtrim(data_content))
    error('data_content cannot be empty.');
end

if ~ismember(regression_mode, {'all', 'central'})
    error('regression_mode must be ''all'' or ''central''.');
end

if ~ismember(slope_plot_style, {'bar_points', 'violin', 'points_median'})
    error(['slope_plot_style must be ''bar_points'', ''violin'', or ' ...
           '''points_median''.']);
end

if ~ismember(delta_plot_style, {'bar_points', 'violin', 'points_median'})
    error(['delta_plot_style must be ''bar_points'', ''violin'', or ' ...
           '''points_median''.']);
end

validatePlotStatisticOptionsLocal(slope_show_mean, 'slope_show_mean');
validatePlotStatisticOptionsLocal(slope_show_median, 'slope_show_median');
validatePlotStatisticOptionsLocal(delta_show_mean, 'delta_show_mean');
validatePlotStatisticOptionsLocal(delta_show_median, 'delta_show_median');

slope_show_mean_current = logical(slope_show_mean.(slope_plot_style));
slope_show_median_current = logical(slope_show_median.(slope_plot_style));
delta_show_mean_current = logical(delta_show_mean.(delta_plot_style));
delta_show_median_current = logical(delta_show_median.(delta_plot_style));

validateattributes(runIdx, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'runIdx');
validateattributes(time_window_index, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, ...
    mfilename, 'time_window_index');
validateattributes(analysis_window, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'}, ...
    mfilename, 'analysis_window');
validateattributes(sliding_step, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'}, ...
    mfilename, 'sliding_step');
if ~(isnumeric(group_colors_base) && isreal(group_colors_base) && ...
        size(group_colors_base, 2) == 3 && ...
        all(isfinite(group_colors_base(:))) && ...
        all(group_colors_base(:) >= 0) && ...
        all(group_colors_base(:) <= 1))
    error('group_colors_base must be an N x 3 numeric matrix in [0, 1].');
end

if ~(isnumeric(session_shade_range) && numel(session_shade_range) == 2 && ...
        all(isfinite(session_shade_range)) && ...
        session_shade_range(1) > 0 && ...
        session_shade_range(2) <= 1 && ...
        session_shade_range(2) >= session_shade_range(1))
    error('session_shade_range must be an increasing two-element vector in (0, 1].');
end

if strcmp(effect_type, 'size')
    validateattributes(pick_contrast, {'numeric'}, ...
        {'scalar', 'integer', 'finite'}, mfilename, 'pick_contrast');
    if ~ismember(pick_contrast, [0, 1, 2])
        error('pick_contrast must be 0, 1, or 2.');
    end
end

validateLogicalScalarLocal(show_session_fit_lines, 'show_session_fit_lines');
validateLogicalScalarLocal(skip_invalid_sessions, 'skip_invalid_sessions');
validateLogicalScalarLocal(save_fig, 'save_fig');
validateLogicalScalarLocal(save_svg, 'save_svg');
validateLogicalScalarLocal(save_png, 'save_png');
validateLogicalScalarLocal(save_mat, 'save_mat');
validateLogicalScalarLocal(close_after_save, 'close_after_save');
validateLogicalScalarLocal(scatter_show_legend, 'scatter_show_legend');
validateLogicalScalarLocal(use_broken_axis_display, ...
    'use_broken_axis_display');
validateLogicalScalarLocal(force_symmetric_broken_axis, ...
    'force_symmetric_broken_axis');
validateLogicalScalarLocal(draw_broken_axis_marks, ...
    'draw_broken_axis_marks');
validateLogicalScalarLocal(use_delta_overflow_display, ...
    'use_delta_overflow_display');

show_session_fit_lines = logical(show_session_fit_lines);
skip_invalid_sessions = logical(skip_invalid_sessions);
save_fig = logical(save_fig);
save_svg = logical(save_svg);
save_png = logical(save_png);
save_any_figure = save_fig || save_svg || save_png;
save_mat = logical(save_mat);
close_after_save = logical(close_after_save);
scatter_show_legend = logical(scatter_show_legend);
use_broken_axis_display = logical(use_broken_axis_display);
force_symmetric_broken_axis = logical(force_symmetric_broken_axis);
draw_broken_axis_marks = logical(draw_broken_axis_marks);
use_delta_overflow_display = logical(use_delta_overflow_display);

validateattributes(png_dpi, {'numeric'}, ...
    {'scalar', 'real', 'integer', 'positive', 'finite'}, ...
    mfilename, 'png_dpi');

validateattributes(delta_overflow_prctile, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 100}, ...
    mfilename, 'delta_overflow_prctile');
validateattributes(delta_overflow_trigger_ratio, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>=', 1}, ...
    mfilename, 'delta_overflow_trigger_ratio');
validateattributes(delta_overflow_gap_frac, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'}, ...
    mfilename, 'delta_overflow_gap_frac');
validateattributes(delta_overflow_top_pad_frac, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'}, ...
    mfilename, 'delta_overflow_top_pad_frac');

if ~(isnumeric(delta_overflow_thresholds) && ...
        numel(delta_overflow_thresholds) == 2 && ...
        all(isfinite(delta_overflow_thresholds) | ...
            isnan(delta_overflow_thresholds)))
    error(['delta_overflow_thresholds must contain two finite values or ' ...
           'NaNs.']);
end
delta_overflow_thresholds = reshape(delta_overflow_thresholds, 1, 2);

group_names = normalizeManualGroupNamesLocal(group_names);
[group_display_names, group_file_tags, group_mapping_file_tag] = ...
    buildGroupLabelsLocal(group_names);
nGroups = numel(group_names);
group_colors = resolveGroupColorsLocal(group_colors_base, nGroups);

if ~(isnumeric(comparison_x_positions) && ...
        numel(comparison_x_positions) == 4 && ...
        all(isfinite(comparison_x_positions)) && ...
        all(diff(comparison_x_positions) > 0))
    error('comparison_x_positions must contain four increasing finite values.');
end
comparison_x_positions = reshape(comparison_x_positions, 1, []);

if ~isempty(slope_y_limits)
    slope_y_limits = validateAxisLimitsLocal(slope_y_limits, 'slope_y_limits');
end

delta_y_limits = validateTwoRowAxisLimitsLocal( ...
    delta_y_limits, 'delta_y_limits');

field_names = { ...
    field_original, ...
    field_all, ...
    field_across, ...
    field_within, ...
    field_feedforward, ...
    field_feedback};

field_names = cellstr(string(field_names(:)))';

field_labels = { ...
    'Original', ...
    'Use all', ...
    'Across', ...
    'Within', ...
    'Feedforward', ...
    'Feedback'};

field_axis_labels = { ...
    'Original data', ...
    'All latents', ...
    'Across latents', ...
    'Within latents', ...
    'Feedforward latents', ...
    'Feedback latents'};

if numel(unique(field_names, 'stable')) ~= numel(field_names)
    error('The six fixed field names must be unique.');
end

comparison_defs = makeComparisonDefinitionsLocal( ...
    field_names, field_labels, field_axis_labels);
effect_metric_label = effectMetricLabelLocal(effect_type, effect_metric);
effect_modulation_label = effectModulationLabelLocal(effect_type);

fprintf('Effect type          : %s\n', effect_type);
fprintf('Effect metric        : %s\n', effect_metric);
fprintf('Data content         : %s\n', data_content);
fprintf('Model mode           : %s\n', model_mode);
fprintf('Run                   : %d\n', runIdx);
fprintf('Time window          : %.12g s, step %.12g s, window %d\n', ...
    analysis_window, sliding_step, time_window_index);
fprintf('Regression mode      : %s\n', regression_mode);
fprintf('Show session fits    : %d\n', show_session_fit_lines);
fprintf('Slope style          : %s\n', slope_plot_style);
fprintf('Delta style          : %s\n', delta_plot_style);
fprintf('Group-area mapping   : %s\n', group_mapping_file_tag);

%% ======================= FIND SESSION FOLDERS ==========================

session_dirs = findCatgtSessionDirsLocal(root_dir);

if isempty(session_dirs)
    error('No catgt_* session folders found one level below root_dir: %s', root_dir);
end

fprintf('Found %d candidate catgt_* session folders.\n', numel(session_dirs));

%% ======================== READ ALL SESSIONS ============================

sessions = struct( ...
    'session_name', {}, ...
    'session_dir', {}, ...
    'source_files', {}, ...
    'groupd', {}, ...
    'stim_tag', {}, ...
    'metric_by_field', {});

skipped = struct('session_dir', {}, 'reason', {});

for si = 1:numel(session_dirs)
    session_dir = session_dirs{si};
    [~, session_name] = fileparts(session_dir);

    fprintf('\n[%d/%d] Reading %s\n', si, numel(session_dirs), session_dir);

    try
        compact_results = cell(1, numel(field_names));
        source_files = cell(1, numel(field_names));

        for f = 1:numel(field_names)
            [compact_results{f}, source_files{f}] = ...
                locateAndLoadEffectResultLocal( ...
                session_dir, effect_type, effect_metric, ...
                data_content, model_mode, runIdx, pick_contrast, ...
                analysis_window, sliding_step, time_window_index, ...
                field_names{f}, group_mapping_file_tag);

            fprintf('  %-12s: %s\n', field_labels{f}, source_files{f});
        end

        validateSessionFieldCompatibilityLocal( ...
            compact_results, field_names, session_name);

        nGroups_this = numel(compact_results{1}.groupd);

        if nGroups_this ~= nGroups
            error(['Session %s contains %d model groups, but group_names ' ...
                   'contains %d labels.'], ...
                  session_name, nGroups_this, nGroups);
        end

        metric_by_field = cell(numel(field_names), nGroups);

        for f = 1:numel(field_names)
            for g = 1:nGroups
                metric_by_field{f, g} = ...
                    compact_results{f}.metric_by_group{g};
            end
        end

        rec = struct();
        rec.session_name = session_name;
        rec.session_dir = session_dir;
        rec.source_files = source_files;
        rec.groupd = compact_results{1}.groupd;
        rec.stim_tag = compact_results{1}.stim_tag;
        rec.metric_by_field = metric_by_field;

        sessions(end + 1) = rec; %#ok<SAGROW>

    catch ME
        reason = ME.message;

        if ~skip_invalid_sessions
            rethrow(ME);
        end

        warning('Skipping %s:\n%s', session_dir, reason);
        skipped(end + 1) = struct( ...
            'session_dir', session_dir, ...
            'reason', reason); %#ok<SAGROW>
    end
end

if isempty(sessions)
    error('No valid sessions contained all six compatible effect result files.');
end

nSessions = numel(sessions);

if isempty(session_display_names)
    session_display_names = arrayfun( ...
        @(k) sprintf('S%d', k), 1:nSessions, 'UniformOutput', false);
else
    session_display_names = cellstr(string(session_display_names(:)))';

    if numel(session_display_names) ~= nSessions
        error(['session_display_names contains %d labels, but %d valid ' ...
               'sessions were loaded.'], ...
              numel(session_display_names), nSessions);
    end
end

[session_colors, session_shade_weights, session_gray_colors] = ...
    buildSessionColorsLocal(group_colors, nSessions, session_shade_range);

fprintf('\nValid sessions and shade order:\n');
for s = 1:nSessions
    fprintf('  %s -> %s (shade %.3f)\n', ...
        session_display_names{s}, sessions(s).session_name, ...
        session_shade_weights(s));
end

%% =================== BUILD SCATTER/REGRESSION DATA =====================

scatter_summary = repmat(struct( ...
    'key', '', ...
    'base_field', '', ...
    'comparison_field', '', ...
    'base_label', '', ...
    'comparison_label', '', ...
    'base_axis_label', '', ...
    'comparison_axis_label', '', ...
    'display_axis_info', struct(), ...
    'session_axis_info', [], ...
    'session_group', []), ...
    1, numel(comparison_defs));

axis_options = struct();
axis_options.use_broken_axis_display = use_broken_axis_display;
axis_options.break_start_prctile = break_start_prctile;
axis_options.broken_axis_trigger_ratio = broken_axis_trigger_ratio;
axis_options.tail_display_frac = tail_display_frac;
axis_options.break_gap_frac = break_gap_frac;
axis_options.min_tail_display_slope = min_tail_display_slope;
axis_options.force_symmetric_broken_axis = force_symmetric_broken_axis;

empty_reg = computeLinearRegressionLocal([], []);
empty_session_group = struct( ...
    'x', [], ...
    'y', [], ...
    'valid_mask', [], ...
    'central_mask_for_valid', [], ...
    'regression_all', empty_reg, ...
    'regression_central', empty_reg);

for p = 1:numel(comparison_defs)
    def = comparison_defs(p);
    all_x = [];
    all_y = [];

    for s = 1:nSessions
        for g = 1:nGroups
            x = sessions(s).metric_by_field{def.base_field_index, g};
            y = sessions(s).metric_by_field{def.comparison_field_index, g};
            valid = isfinite(x) & isfinite(y);
            all_x = [all_x; x(valid)]; %#ok<AGROW>
            all_y = [all_y; y(valid)]; %#ok<AGROW>
        end
    end

    display_axis_info = computeSymmetricBrokenAxisInfoLocal( ...
        all_x, all_y, axis_options);

    session_axis_info = repmat(display_axis_info, nSessions, 1);
    session_group = repmat(empty_session_group, nSessions, nGroups);

    for s = 1:nSessions
        session_x = [];
        session_y = [];

        for g = 1:nGroups
            x = sessions(s).metric_by_field{def.base_field_index, g};
            y = sessions(s).metric_by_field{def.comparison_field_index, g};
            valid = isfinite(x) & isfinite(y);
            session_x = [session_x; x(valid)]; %#ok<AGROW>
            session_y = [session_y; y(valid)]; %#ok<AGROW>
        end

        session_axis_info(s) = computeSymmetricBrokenAxisInfoLocal( ...
            session_x, session_y, axis_options);

        for g = 1:nGroups
            x_raw = sessions(s).metric_by_field{def.base_field_index, g};
            y_raw = sessions(s).metric_by_field{def.comparison_field_index, g};
            valid = isfinite(x_raw) & isfinite(y_raw);
            x = x_raw(valid);
            y = y_raw(valid);

            central = computeCentralMaskLocal( ...
                x, y, session_axis_info(s));

            session_group(s, g).x = x;
            session_group(s, g).y = y;
            session_group(s, g).valid_mask = valid;
            session_group(s, g).central_mask_for_valid = central;
            session_group(s, g).regression_all = ...
                computeLinearRegressionLocal(x, y);
            session_group(s, g).regression_central = ...
                computeLinearRegressionLocal(x(central), y(central));
        end
    end

    scatter_summary(p).key = def.key;
    scatter_summary(p).base_field = def.base_field;
    scatter_summary(p).comparison_field = def.comparison_field;
    scatter_summary(p).base_label = def.base_label;
    scatter_summary(p).comparison_label = def.comparison_label;
    scatter_summary(p).base_axis_label = def.base_axis_label;
    scatter_summary(p).comparison_axis_label = def.comparison_axis_label;
    scatter_summary(p).display_axis_info = display_axis_info;
    scatter_summary(p).session_axis_info = session_axis_info;
    scatter_summary(p).session_group = session_group;
end

%% ========================== BUILD SLOPES ===============================

nComponentComparisons = 4;
slope_values = nan(nSessions, nComponentComparisons, nGroups);
slope_regression = repmat( ...
    empty_reg, nSessions, nComponentComparisons, nGroups);

for c = 1:nComponentComparisons
    scatter_panel_idx = c + 1;
    this_panel = scatter_summary(scatter_panel_idx);

    for s = 1:nSessions
        for g = 1:nGroups
            if strcmp(regression_mode, 'central')
                reg = this_panel.session_group(s, g).regression_central;
            else
                reg = this_panel.session_group(s, g).regression_all;
            end

            slope_regression(s, c, g) = reg;

            if reg.is_valid
                slope_values(s, c, g) = reg.slope;
            end
        end
    end
end

slope_stats = repmat(computeBasicStatsLocal([]), ...
    nComponentComparisons, nGroups);

for c = 1:nComponentComparisons
    for g = 1:nGroups
        slope_stats(c, g) = computeBasicStatsLocal( ...
            slope_values(:, c, g));
    end
end

%% ========================== BUILD DELTAS ===============================

% Dimensions: comparison x group x base-sign x session.
% base-sign index 1 = base >= 0; index 2 = base < 0.
delta_by_session = cell(nComponentComparisons, nGroups, 2, nSessions);
delta_pooled = cell(nComponentComparisons, nGroups, 2);
delta_stats = repmat(computeBasicStatsLocal([]), ...
    nComponentComparisons, nGroups, 2);

for c = 1:nComponentComparisons
    scatter_panel_idx = c + 1;

    for g = 1:nGroups
        for sign_idx = 1:2
            pooled = [];

            for s = 1:nSessions
                x = scatter_summary(scatter_panel_idx).session_group(s, g).x;
                y = scatter_summary(scatter_panel_idx).session_group(s, g).y;

                if sign_idx == 1
                    selected = x >= 0;
                else
                    selected = x < 0;
                end

                delta = y(selected) - x(selected);
                delta = delta(isfinite(delta));
                delta_by_session{c, g, sign_idx, s} = delta(:);
                pooled = [pooled; delta(:)]; %#ok<AGROW>
            end

            delta_pooled{c, g, sign_idx} = pooled;
            delta_stats(c, g, sign_idx) = computeBasicStatsLocal(pooled);
        end
    end
end

delta_overflow_options = struct();
delta_overflow_options.enabled = use_delta_overflow_display;
delta_overflow_options.percentile = delta_overflow_prctile;
delta_overflow_options.trigger_ratio = delta_overflow_trigger_ratio;
delta_overflow_options.manual_thresholds = delta_overflow_thresholds;
delta_overflow_options.gap_frac = delta_overflow_gap_frac;
delta_overflow_options.top_pad_frac = delta_overflow_top_pad_frac;
delta_overflow_options.y_limits = delta_y_limits;

delta_overflow_info = computeDeltaOverflowInfoLocal( ...
    delta_pooled, delta_overflow_options);

%% ======================= OUTPUT NAMES AND SUMMARY ======================

effect_tag = effectTypeFileTagLocal(effect_type);
metric_tag = effectMetricFileTagLocal(effect_metric);
model_tag = modelModeFileTagLocal(model_mode);
selection_tag = effectSelectionTagLocal(effect_type, pick_contrast);
time_tag = sprintf('wn%s_st%s_w%02d', ...
    formatSecondsForFilenameLocal(analysis_window), ...
    formatSecondsForFilenameLocal(sliding_step), ...
    time_window_index);

if show_session_fit_lines
    fit_tag = ['fit_', regression_mode];
else
    fit_tag = 'nofit';
end

if slope_show_mean_current
    slope_mean_tag = 'mean';
else
    slope_mean_tag = 'nomean';
end

if slope_show_median_current
    slope_median_tag = 'median';
else
    slope_median_tag = 'nomedian';
end

if delta_show_mean_current
    delta_mean_tag = 'mean';
else
    delta_mean_tag = 'nomean';
end

if delta_show_median_current
    delta_median_tag = 'median';
else
    delta_median_tag = 'nomedian';
end


delta_overflow_tag = deltaOverflowFileTagLocal( ...
    use_delta_overflow_display, delta_overflow_prctile, ...
    delta_overflow_thresholds);

out_base = sprintf( ...
    ['sum_sessions_effect_%s_%s_%s_run%03d_%s_%s_%s_%s_%s_' ...
     'slope_%s_%s_%s_delta_%s_%s_%s_%s'], ...
    effect_tag, sanitizeFilenameLocal(data_content), model_tag, runIdx, ...
    group_mapping_file_tag, metric_tag, selection_tag, time_tag, fit_tag, ...
    slope_plot_style, slope_mean_tag, slope_median_tag, ...
    delta_plot_style, delta_mean_tag, delta_median_tag, ...
    delta_overflow_tag);

figure_suffixes = { ...
    '01_original_vs_all', ...
    '02_all_vs_across', ...
    '03_all_vs_within', ...
    '04_all_vs_feedforward', ...
    '05_all_vs_feedback', ...
    '06_slope', ...
    '07_delta'};

figure_files = repmat(struct('fig', '', 'svg', '', 'png', ''), 1, 7);

for k = 1:7
    figure_files(k).fig = fullfile( ...
        root_dir, sprintf('%s_%s.fig', out_base, figure_suffixes{k}));
    figure_files(k).svg = fullfile( ...
        root_dir, sprintf('%s_%s.svg', out_base, figure_suffixes{k}));
    figure_files(k).png = fullfile( ...
        root_dir, sprintf('%s_%s.png', out_base, figure_suffixes{k}));
end

mat_file = fullfile(root_dir, sprintf('%s_summary.mat', out_base));

EffectSummary = struct();
EffectSummary.meta = struct();
EffectSummary.meta.program = mfilename;
EffectSummary.meta.revision = '2026-09-08a';
EffectSummary.meta.root_dir = root_dir;
EffectSummary.meta.effect_type = effect_type;
EffectSummary.meta.effect_metric = effect_metric;
EffectSummary.meta.effect_metric_label = effect_metric_label;
EffectSummary.meta.effect_modulation_label = effect_modulation_label;
EffectSummary.meta.data_content = data_content;
EffectSummary.meta.model_mode = model_mode;
EffectSummary.meta.runIdx = runIdx;
EffectSummary.meta.pick_contrast = pick_contrast;
EffectSummary.meta.analysis_window = analysis_window;
EffectSummary.meta.sliding_step = sliding_step;
EffectSummary.meta.time_window_index = time_window_index;
EffectSummary.meta.regression_mode = regression_mode;
EffectSummary.meta.show_session_fit_lines = show_session_fit_lines;
EffectSummary.meta.use_broken_axis_display = use_broken_axis_display;
EffectSummary.meta.break_start_prctile = break_start_prctile;
EffectSummary.meta.broken_axis_trigger_ratio = broken_axis_trigger_ratio;
EffectSummary.meta.tail_display_frac = tail_display_frac;
EffectSummary.meta.break_gap_frac = break_gap_frac;
EffectSummary.meta.min_tail_display_slope = min_tail_display_slope;
EffectSummary.meta.force_symmetric_broken_axis = force_symmetric_broken_axis;
EffectSummary.meta.slope_plot_style = slope_plot_style;
EffectSummary.meta.slope_show_mean = slope_show_mean;
EffectSummary.meta.slope_show_median = slope_show_median;
EffectSummary.meta.slope_show_mean_current = slope_show_mean_current;
EffectSummary.meta.slope_show_median_current = slope_show_median_current;
EffectSummary.meta.delta_plot_style = delta_plot_style;
EffectSummary.meta.delta_show_mean = delta_show_mean;
EffectSummary.meta.delta_show_median = delta_show_median;
EffectSummary.meta.delta_show_mean_current = delta_show_mean_current;
EffectSummary.meta.delta_show_median_current = delta_show_median_current;
EffectSummary.meta.use_delta_overflow_display = ...
    use_delta_overflow_display;
EffectSummary.meta.delta_overflow_prctile = delta_overflow_prctile;
EffectSummary.meta.delta_overflow_trigger_ratio = ...
    delta_overflow_trigger_ratio;
EffectSummary.meta.delta_overflow_thresholds = ...
    delta_overflow_thresholds;
EffectSummary.meta.delta_overflow_gap_frac = delta_overflow_gap_frac;
EffectSummary.meta.delta_overflow_top_pad_frac = ...
    delta_overflow_top_pad_frac;
EffectSummary.meta.delta_overflow_display_rule = ...
    ['Values above each sign row threshold are shown at one common top ' ...
     'overflow level. Raw values and raw statistics are unchanged.'];
EffectSummary.meta.delta_mean_rule = ...
    ['Pool every finite neuron delta from every session directly, then ' ...
     'calculate the arithmetic mean; sessions are not averaged first.'];
EffectSummary.meta.delta_definition = ...
    'comparison effect metric minus Use-all effect metric';
EffectSummary.meta.base_sign_rule = ...
    'positive base uses Use-all >= 0; negative base uses Use-all < 0';
EffectSummary.meta.num_candidate_sessions = numel(session_dirs);
EffectSummary.meta.num_valid_sessions = nSessions;
EffectSummary.meta.num_skipped_sessions = numel(skipped);
EffectSummary.meta.session_display_names = session_display_names;
EffectSummary.meta.session_shade_weights = session_shade_weights;
EffectSummary.meta.num_groups = nGroups;
EffectSummary.meta.group_names = group_names;
EffectSummary.meta.group_display_names = group_display_names;
EffectSummary.meta.group_file_tags = group_file_tags;
EffectSummary.meta.group_mapping_file_tag = group_mapping_file_tag;
EffectSummary.meta.group_label_source = 'manual group_names parameter';
EffectSummary.meta.group_colors_base = group_colors_base;
EffectSummary.meta.group_colors = group_colors;
EffectSummary.meta.session_colors = session_colors;
EffectSummary.meta.save_fig = save_fig;
EffectSummary.meta.save_svg = save_svg;
EffectSummary.meta.save_png = save_png;
EffectSummary.meta.png_dpi = png_dpi;
EffectSummary.meta.save_any_figure = save_any_figure;
EffectSummary.meta.save_mat = save_mat;

EffectSummary.field_names = field_names;
EffectSummary.field_labels = field_labels;
EffectSummary.field_axis_labels = field_axis_labels;
EffectSummary.comparison_definitions = comparison_defs;
EffectSummary.sessions = sessions;
EffectSummary.skipped = skipped;
EffectSummary.scatter = scatter_summary;
EffectSummary.slope.values = slope_values;
EffectSummary.slope.regression = slope_regression;
EffectSummary.slope.stats = slope_stats;
EffectSummary.slope.comparison_labels = field_labels(3:6);
EffectSummary.slope.comparison_x_positions = comparison_x_positions;
EffectSummary.delta.values_by_session = delta_by_session;
EffectSummary.delta.pooled_values = delta_pooled;
EffectSummary.delta.stats = delta_stats;
EffectSummary.delta.overflow_info = delta_overflow_info;
EffectSummary.delta.comparison_labels = field_labels(3:6);
EffectSummary.delta.base_sign_labels = {'base_ge_0', 'base_lt_0'};
EffectSummary.output_files.mat = mat_file;
EffectSummary.output_files.figures = figure_files;

if save_mat
    save(mat_file, 'EffectSummary', '-v7.3');
    fprintf('\nSaved summary mat:\n  %s\n', mat_file);
else
    fprintf('\nMAT saving is disabled.\n');
end

%% ============================== PLOT ===================================

plot_config = struct();
plot_config.figure_visible = figure_visible;
plot_config.font_name = font_name;
plot_config.font_size = font_size;
plot_config.axis_line_width = axis_line_width;
plot_config.group_colors = group_colors;
plot_config.group_display_names = group_display_names;
plot_config.session_colors = session_colors;
plot_config.session_gray_colors = session_gray_colors;
plot_config.session_display_names = session_display_names;
plot_config.effect_metric_label = effect_metric_label;
plot_config.effect_modulation_label = effect_modulation_label;
plot_config.regression_mode = regression_mode;
plot_config.show_session_fit_lines = show_session_fit_lines;
plot_config.draw_broken_axis_marks = draw_broken_axis_marks;
plot_config.scatter_point_size = scatter_point_size;
plot_config.scatter_point_alpha = scatter_point_alpha;
plot_config.fit_line_width = fit_line_width;
plot_config.reference_line_width = reference_line_width;
plot_config.zero_line_width = zero_line_width;
plot_config.scatter_figure_size_inches = scatter_figure_size_inches;
plot_config.scatter_show_legend = scatter_show_legend;

figure_handles = gobjects(1, 7);

for p = 1:5
    allow_fit = p > 1 && show_session_fit_lines;
    figure_handles(p) = plotScatterSummaryLocal( ...
        scatter_summary(p), allow_fit, plot_config);

    if save_any_figure
        saveFigureFilesLocal(figure_handles(p), figure_files(p), ...
            save_fig, save_svg, save_png, png_dpi);
    end
end

slope_config = plot_config;
slope_config.figure_size_inches = slope_figure_size_inches;
slope_config.plot_style = slope_plot_style;
slope_config.show_mean = slope_show_mean_current;
slope_config.show_median = slope_show_median_current;
slope_config.y_limits = slope_y_limits;
slope_config.bar_width = slope_bar_width;
slope_config.bar_alpha = slope_bar_alpha;
slope_config.point_size = slope_point_size;
slope_config.point_alpha = slope_point_alpha;
slope_config.point_jitter_width = slope_point_jitter_width;
slope_config.violin_width = slope_violin_width;
slope_config.violin_alpha = slope_violin_alpha;
slope_config.median_line_half_width = slope_median_line_half_width;
slope_config.mean_marker_size = slope_mean_marker_size;
slope_config.group_offset = slope_group_offset;
slope_config.x_positions = comparison_x_positions;
slope_config.comparison_labels = field_labels(3:6);

figure_handles(6) = plotSlopeSummaryLocal(slope_values, slope_config);

if save_any_figure
    saveFigureFilesLocal(figure_handles(6), figure_files(6), ...
        save_fig, save_svg, save_png, png_dpi);
end

delta_config = plot_config;
delta_config.figure_size_inches = delta_figure_size_inches;
delta_config.plot_style = delta_plot_style;
delta_config.show_mean = delta_show_mean_current;
delta_config.show_median = delta_show_median_current;
delta_config.overflow_info = delta_overflow_info;
delta_config.bar_width = delta_bar_width;
delta_config.bar_alpha = delta_bar_alpha;
delta_config.point_size = delta_point_size;
delta_config.point_alpha = delta_point_alpha;
delta_config.point_jitter_width = delta_point_jitter_width;
delta_config.violin_width = delta_violin_width;
delta_config.violin_alpha = delta_violin_alpha;
delta_config.median_line_half_width = delta_median_line_half_width;
delta_config.mean_marker_size = delta_mean_marker_size;
delta_config.group_offset = delta_group_offset;
delta_config.x_positions = comparison_x_positions;
delta_config.comparison_labels = field_labels(3:6);
delta_config.effect_type = effect_type;

figure_handles(7) = plotDeltaSummaryLocal( ...
    delta_by_session, delta_pooled, delta_stats, delta_config);

if save_any_figure
    saveFigureFilesLocal(figure_handles(7), figure_files(7), ...
        save_fig, save_svg, save_png, png_dpi);
end

if save_any_figure
    saved_format_names = {};
    if save_fig
        saved_format_names{end + 1} = 'FIG'; %#ok<SAGROW>
    end
    if save_svg
        saved_format_names{end + 1} = 'SVG'; %#ok<SAGROW>
    end
    if save_png
        saved_format_names{end + 1} = sprintf( ...
            'PNG (%d dpi)', png_dpi); %#ok<SAGROW>
    end

    fprintf('\nSaved seven figure sets [%s] in:\n  %s\n', ...
        strjoin(saved_format_names, ', '), root_dir);
else
    fprintf(['\nFigure saving is disabled because save_fig, save_svg, ' ...
             'and save_png are all false.\n']);
end

if ~isempty(skipped)
    fprintf('Skipped %d session folder(s); inspect warnings', numel(skipped));

    if save_mat
        fprintf(' or EffectSummary.skipped.\n');
    else
        fprintf('.\n');
    end
end

if close_after_save
    close(figure_handles(isgraphics(figure_handles)));
end

fprintf('Done.\n');

%% ======================================================================
% Local functions
% =======================================================================

function effect_type = normalizeEffectTypeLocal(effect_type)

    effect_type = lower(strtrim(char(effect_type)));

    switch regexprep(effect_type, '[^a-z]', '')
        case {'size', 'sizeeffect'}
            effect_type = 'size';

        case {'contrast', 'contrasteffect'}
            effect_type = 'contrast';

        otherwise
            error('effect_type must be ''size'' or ''contrast''.');
    end
end

function model_mode = normalizeModelModeLocal(model_mode)

    model_mode = lower(strtrim(char(model_mode)));
    compact = regexprep(model_mode, '[^a-z0-9]', '');

    switch compact
        case {'allconditionmodel', 'allconditionsmodel', 'allcondition'}
            model_mode = 'all_condition_model';

        case {'conditionspecificmodels', 'conditionspecificmodel', ...
              'conditionmodels', 'conditionmode'}
            model_mode = 'condition_specific_models';

        otherwise
            error(['Unknown model_mode: %s. Use ''all_condition_model'' ' ...
                   'or ''condition_specific_models''.'], model_mode);
    end
end

function metric = resolveEffectMetricLocal(effect_type, metric)

    metric = char(string(metric));
    metric = strtrim(metric);

    if strcmpi(metric, 'auto')
        if strcmp(effect_type, 'size')
            metric = 'S_norm_diff';
        else
            metric = 'CMI';
        end
    end

    if strcmp(effect_type, 'size')
        valid = {'classic_SI', 'delta_SL', 'S_norm_diff'};
    else
        valid = {'delta_HL', 'CMI'};
    end

    match = find(strcmpi(metric, valid), 1);

    if isempty(match)
        error('Metric %s is not valid for effect_type %s.', metric, effect_type);
    end

    metric = valid{match};
end

function validatePlotStatisticOptionsLocal(options, variable_name)

    required_fields = {'bar_points', 'violin', 'points_median'};

    if ~(isstruct(options) && isscalar(options))
        error('%s must be a scalar struct.', variable_name);
    end

    for k = 1:numel(required_fields)
        field_name = required_fields{k};

        if ~isfield(options, field_name)
            error('%s is missing field %s.', variable_name, field_name);
        end

        validateLogicalScalarLocal( ...
            options.(field_name), sprintf('%s.%s', variable_name, field_name));
    end
end

function validateLogicalScalarLocal(value, variable_name)

    valid = (islogical(value) || isnumeric(value)) && ...
        isscalar(value) && isfinite(double(value)) && ...
        ismember(double(value), [0, 1]);

    if ~valid
        error('%s must be true or false.', variable_name);
    end
end

function limits = validateAxisLimitsLocal(limits, variable_name)

    if ~(isnumeric(limits) && numel(limits) == 2 && ...
            all(isfinite(limits)) && limits(2) > limits(1))
        error('%s must be a finite increasing two-element vector.', ...
            variable_name);
    end

    limits = reshape(limits, 1, 2);
end

function limits = validateTwoRowAxisLimitsLocal(limits, variable_name)

    if isempty(limits)
        limits = [];
        return;
    end

    if ~(isnumeric(limits) && all(isfinite(limits(:))))
        error('%s must be empty or contain finite numeric limits.', ...
            variable_name);
    end

    if isvector(limits) && numel(limits) == 2
        one_row = validateAxisLimitsLocal(limits, variable_name);
        limits = repmat(one_row, 2, 1);
        return;
    end

    if ~isequal(size(limits), [2, 2])
        error('%s must be empty, 1 x 2, or 2 x 2.', variable_name);
    end

    for row_idx = 1:2
        if limits(row_idx, 2) <= limits(row_idx, 1)
            error('Each row of %s must be increasing.', variable_name);
        end
    end
end

function defs = makeComparisonDefinitionsLocal( ...
        field_names, field_labels, field_axis_labels)

    base_indices = [1, 2, 2, 2, 2];
    comparison_indices = [2, 3, 4, 5, 6];
    keys = { ...
        'original_vs_all', ...
        'all_vs_across', ...
        'all_vs_within', ...
        'all_vs_feedforward', ...
        'all_vs_feedback'};

    defs = repmat(struct( ...
        'key', '', ...
        'base_field_index', [], ...
        'comparison_field_index', [], ...
        'base_field', '', ...
        'comparison_field', '', ...
        'base_label', '', ...
        'comparison_label', '', ...
        'base_axis_label', '', ...
        'comparison_axis_label', ''), ...
        1, numel(keys));

    for k = 1:numel(keys)
        b = base_indices(k);
        c = comparison_indices(k);
        defs(k).key = keys{k};
        defs(k).base_field_index = b;
        defs(k).comparison_field_index = c;
        defs(k).base_field = field_names{b};
        defs(k).comparison_field = field_names{c};
        defs(k).base_label = field_labels{b};
        defs(k).comparison_label = field_labels{c};
        defs(k).base_axis_label = field_axis_labels{b};
        defs(k).comparison_axis_label = field_axis_labels{c};
    end
end

function label = effectModulationLabelLocal(effect_type)

    if strcmp(effect_type, 'size')
        label = 'Size modulation';
    else
        label = 'Contrast modulation';
    end
end

function label = effectMetricLabelLocal(effect_type, metric)

    if strcmp(effect_type, 'size')
        switch metric
            case 'classic_SI'
                label = 'Size index';

            case 'delta_SL'
                label = 'S - L response';

            case 'S_norm_diff'
                label = 'Normalized size effect';
        end
    else
        switch metric
            case 'delta_HL'
                label = 'H - L response';

            case 'CMI'
                label = 'Contrast modulation index';
        end
    end
end

function session_dirs = findCatgtSessionDirsLocal(root_dir)

    session_dirs = {};
    [~, root_name] = fileparts(root_dir);

    if startsWith(lower(root_name), 'catgt_')
        session_dirs = {root_dir};
    end

    listing = dir(fullfile(root_dir, '*', 'catgt_*'));

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

function [compact_result, selected_file] = ...
        locateAndLoadEffectResultLocal( ...
        session_dir, effect_type, effect_metric, ...
        data_content, model_mode, runIdx, pick_contrast, ...
        analysis_window, sliding_step, time_window_index, analysis_field, ...
        group_mapping_file_tag)

    safe_data = sanitizeFilenameLocal(data_content);
    safe_field = sanitizeFilenameLocal(analysis_field);

    base_name = sprintf('%s_%s_%s_effect_%s', ...
        safe_data, model_mode, effect_type, safe_field);

    if strcmp(effect_type, 'size')
        base_name = [base_name, sizeContrastFileSuffixLocal(pick_contrast)];
    end

    area_tagged_base_name = sprintf( ...
        '%s_%s', base_name, group_mapping_file_tag);
    time_parameter_tag = sprintf('_wn%s_st%s', ...
        formatSecondsForFilenameLocal(analysis_window), ...
        formatSecondsForFilenameLocal(sliding_step));

    tagged_name = sprintf('%s%s_w%02d.mat', ...
        area_tagged_base_name, time_parameter_tag, time_window_index);
    names_to_find = {tagged_name};

    if time_window_index == 1
        % A full-trial single-window output has no time suffix, but it always
        % carries the manually specified group-area mapping tag.
        names_to_find{end + 1} = sprintf( ...
            '%s.mat', area_tagged_base_name); %#ok<AGROW>
    end

    if strcmp(model_mode, 'all_condition_model')
        preferred_dir = fullfile( ...
            session_dir, ['FA_Dlag_', data_content], ...
            'mat_results', sprintf('run%03d', runIdx));
    else
        preferred_dir = session_dir;
    end

    candidate_files = {};

    % Preferred exact locations are checked first.
    for n = 1:numel(names_to_find)
        f = fullfile(preferred_dir, names_to_find{n});

        if isfile(f)
            candidate_files{end + 1} = f; %#ok<AGROW>
        end
    end

    % Then search recursively within this one session. Exact filenames keep
    % comparison-summary mats and other fields out of the candidate set.
    for n = 1:numel(names_to_find)
        listing = dir(fullfile(session_dir, '**', names_to_find{n}));

        for k = 1:numel(listing)
            if ~listing(k).isdir
                candidate_files{end + 1} = fullfile( ...
                    listing(k).folder, listing(k).name); %#ok<AGROW>
            end
        end
    end

    candidate_files = unique(candidate_files, 'stable');

    if isempty(candidate_files)
        error(['No %s effect result found for field %s and group mapping %s. ' ...
               'Expected one of these filenames: %s within %s.'], ...
              effect_type, analysis_field, group_mapping_file_tag, ...
              strjoin(names_to_find, ', '), session_dir);
    end

    failure_messages = cell(1, numel(candidate_files));
    valid_results = cell(1, numel(candidate_files));
    valid_mask = false(1, numel(candidate_files));

    for k = 1:numel(candidate_files)
        try
            valid_results{k} = loadAndValidateEffectResultLocal( ...
                candidate_files{k}, effect_type, effect_metric, ...
                data_content, model_mode, runIdx, pick_contrast, ...
                analysis_window, sliding_step, time_window_index, ...
                analysis_field);
            valid_mask(k) = true;

        catch ME
            failure_messages{k} = sprintf('%s -> %s', ...
                candidate_files{k}, ME.message);
        end
    end

    valid_idx = find(valid_mask);

    if numel(valid_idx) == 1
        compact_result = valid_results{valid_idx};
        selected_file = candidate_files{valid_idx};
        return;
    end

    if numel(valid_idx) > 1
        valid_files = candidate_files(valid_idx);
        error(['More than one metadata-compatible effect result was found ' ...
               'for field %s and group mapping %s:\n%s'], ...
              analysis_field, group_mapping_file_tag, ...
              strjoin(valid_files, '\n'));
    end

    failure_messages = failure_messages(~cellfun(@isempty, failure_messages));

    error(['Candidate file(s) were found for field %s, but none matched ' ...
           'the requested metadata:\n%s'], ...
          analysis_field, strjoin(failure_messages, '\n'));
end

function compact = loadAndValidateEffectResultLocal( ...
        input_file, effect_type, effect_metric, ...
        data_content, model_mode, runIdx, pick_contrast, ...
        analysis_window, sliding_step, time_window_index, analysis_field)

    if strcmp(effect_type, 'size')
        variable_name = 'size_effect_result';
    else
        variable_name = 'contrast_effect_result';
    end

    S = load(input_file, variable_name);

    if ~isfield(S, variable_name)
        error('%s does not contain %s.', input_file, variable_name);
    end

    result = S.(variable_name);

    if ~isstruct(result) || ~isscalar(result)
        error('%s in %s must be a scalar struct.', variable_name, input_file);
    end

    required_meta = { ...
        'data_content', 'model_mode', 'runIdx', 'analysis_field', ...
        'analysis_window', 'sliding_step', 'time_window_index', ...
        'groupd', effect_metric};

    for k = 1:numel(required_meta)
        if ~isfield(result, required_meta{k})
            error('%s is missing field %s.', variable_name, required_meta{k});
        end
    end

    if ~strcmp(char(string(result.data_content)), data_content)
        error('data_content is %s, expected %s.', ...
            char(string(result.data_content)), data_content);
    end

    result_model_mode = normalizeModelModeLocal(result.model_mode);

    if ~strcmp(result_model_mode, model_mode)
        error('model_mode is %s, expected %s.', ...
            result_model_mode, model_mode);
    end

    if double(result.runIdx) ~= double(runIdx)
        error('runIdx is %g, expected %g.', double(result.runIdx), double(runIdx));
    end

    if ~strcmp(char(string(result.analysis_field)), analysis_field)
        error('analysis_field is %s, expected %s.', ...
            char(string(result.analysis_field)), analysis_field);
    end

    assertNearlyEqualLocal( ...
        result.analysis_window, analysis_window, 'analysis_window');
    assertNearlyEqualLocal( ...
        result.sliding_step, sliding_step, 'sliding_step');

    if double(result.time_window_index) ~= double(time_window_index)
        error('time_window_index is %g, expected %g.', ...
            double(result.time_window_index), double(time_window_index));
    end

    if strcmp(effect_type, 'size')
        if ~isfield(result, 'pick_contrast')
            error('size_effect_result is missing pick_contrast.');
        end

        if double(result.pick_contrast) ~= double(pick_contrast)
            error('pick_contrast is %g, expected %g.', ...
                double(result.pick_contrast), double(pick_contrast));
        end
    end

    groupd = double(result.groupd(:)');

    if isempty(groupd) || any(~isfinite(groupd)) || ...
            any(groupd <= 0) || any(mod(groupd, 1) ~= 0)
        error('groupd must contain finite positive integer group sizes.');
    end

    values = result.(effect_metric);

    if ~(isnumeric(values) && isreal(values) && isvector(values))
        error('%s must be a real numeric vector.', effect_metric);
    end

    values = double(values(:));

    if numel(values) ~= sum(groupd)
        error('%s has %d values, but sum(groupd) is %d.', ...
            effect_metric, numel(values), sum(groupd));
    end

    row_start = 1;
    nGroups = numel(groupd);
    metric_by_group = cell(1, nGroups);

    for g = 1:nGroups
        row_end = row_start + groupd(g) - 1;
        metric_by_group{g} = values(row_start:row_end);
        row_start = row_end + 1;
    end

    compact = struct();
    compact.input_file = input_file;
    compact.analysis_field = analysis_field;
    compact.effect_metric = effect_metric;
    compact.groupd = groupd;
    compact.metric_by_group = metric_by_group;
    compact.analysis_window = double(result.analysis_window);
    compact.sliding_step = double(result.sliding_step);
    compact.time_window_index = double(result.time_window_index);

    if isfield(result, 'stim_tag')
        compact.stim_tag = char(string(result.stim_tag));
    else
        compact.stim_tag = '';
    end
end

function assertNearlyEqualLocal(actual, expected, variable_name)

    if ~(isnumeric(actual) && isscalar(actual) && isfinite(actual))
        error('%s must be a finite numeric scalar.', variable_name);
    end

    tolerance = 1e-9 * max([1, abs(double(actual)), abs(double(expected))]);

    if abs(double(actual) - double(expected)) > tolerance
        error('%s is %.12g, expected %.12g.', ...
            variable_name, double(actual), double(expected));
    end
end

function validateSessionFieldCompatibilityLocal( ...
        compact_results, field_names, session_name)

    reference = compact_results{1};

    for f = 2:numel(compact_results)
        current = compact_results{f};

        if ~isequal(current.groupd, reference.groupd)
            error(['Session %s has a groupd mismatch between %s and %s: ' ...
                   '%s vs %s.'], ...
                  session_name, field_names{1}, field_names{f}, ...
                  mat2str(reference.groupd), mat2str(current.groupd));
        end

        if ~isempty(reference.stim_tag) && ~isempty(current.stim_tag) && ...
                ~strcmp(reference.stim_tag, current.stim_tag)
            error(['Session %s has a stim_tag mismatch between %s and %s: ' ...
                   '%s vs %s.'], ...
                  session_name, field_names{1}, field_names{f}, ...
                  reference.stim_tag, current.stim_tag);
        end
    end
end

function suffix = sizeContrastFileSuffixLocal(pick_contrast)

    switch pick_contrast
        case 0
            suffix = '';

        case 1
            suffix = '_low_contrast';

        case 2
            suffix = '_high_contrast';

        otherwise
            error('pick_contrast must be 0, 1, or 2.');
    end
end

function s = formatSecondsForFilenameLocal(x)

    s = sprintf('%.12g', x);
    s = strrep(s, '.', '_');
    s = strrep(s, '-', 'm');
    s = regexprep(s, '[^a-zA-Z0-9_]', '_');
    s = regexprep(s, '_+', '_');
end

function name = sanitizeFilenameLocal(name)

    name = char(string(name));
    name = regexprep(name, '[^a-zA-Z0-9_\-]', '_');
    name = regexprep(name, '_+', '_');
    name = strtrim(name);

    if isempty(name)
        name = 'unnamed';
    end
end

function group_names = normalizeManualGroupNamesLocal(group_names)

    if isempty(group_names)
        error(['group_names must be specified manually with one area label ' ...
               'for each model group.']);
    end

    if ischar(group_names)
        group_names = {group_names};
    elseif isstring(group_names)
        group_names = cellstr(group_names(:))';
    elseif iscell(group_names)
        group_names = reshape(group_names, 1, []);
    else
        error('group_names must be a char, string array, or cell array.');
    end

    for g = 1:numel(group_names)
        if isstring(group_names{g}) && isscalar(group_names{g})
            group_names{g} = char(group_names{g});
        end

        if ~ischar(group_names{g})
            error('group_names{%d} must be a char or scalar string.', g);
        end

        group_names{g} = strtrim(group_names{g});

        if isempty(group_names{g})
            error('group_names{%d} is empty.', g);
        end
    end
end

function [group_display_names, group_file_tags, group_mapping_file_tag] = ...
        buildGroupLabelsLocal(group_names)

    nGroups = numel(group_names);
    group_display_names = cell(1, nGroups);
    group_file_tags = cell(1, nGroups);

    for g = 1:nGroups
        group_display_names{g} = sprintf( ...
            'Group %d: %s', g, group_names{g});
        group_file_tags{g} = sprintf( ...
            'G%02d_%s', g, sanitizeFilenameLocal(group_names{g}));
    end

    group_mapping_file_tag = strjoin(group_file_tags, '_');
end

function group_colors = resolveGroupColorsLocal(group_colors_base, nGroups)

    automatic_colors = lines(max(nGroups, 1));
    group_colors = automatic_colors(1:nGroups, :);
    nProvided = min(size(group_colors_base, 1), nGroups);

    if nProvided > 0
        group_colors(1:nProvided, :) = group_colors_base(1:nProvided, :);
    end
end

function offsets = centeredGroupOffsetsLocal(nGroups, two_group_half_offset)

    if nGroups <= 1
        offsets = 0;
        return;
    end

    spacing = 2 .* two_group_half_offset;
    offsets = ((1:nGroups) - (nGroups + 1) ./ 2) .* spacing;
end

function tag = effectTypeFileTagLocal(effect_type)

    tag = sanitizeFilenameLocal(effect_type);
end

function tag = modelModeFileTagLocal(model_mode)

    switch model_mode
        case 'all_condition_model'
            tag = 'all';

        case 'condition_specific_models'
            tag = 'csm';
    end
end

function tag = effectMetricFileTagLocal(metric)

    switch metric
        case 'classic_SI'
            tag = 'SI';

        case 'delta_SL'
            tag = 'dSL';

        case 'S_norm_diff'
            tag = 'SnD';

        case 'delta_HL'
            tag = 'dHL';

        case 'CMI'
            tag = 'CMI';

        otherwise
            tag = sanitizeFilenameLocal(metric);
    end
end

function tag = effectSelectionTagLocal(effect_type, pick_contrast)

    if strcmp(effect_type, 'contrast')
        tag = 'relative_contrast';
        return;
    end

    switch pick_contrast
        case 0
            tag = 'all_contrasts';

        case 1
            tag = 'low_contrast';

        case 2
            tag = 'high_contrast';
    end
end

function tag = deltaOverflowFileTagLocal(enabled, percentile, thresholds)

    if ~enabled
        tag = 'dovf_off';
        return;
    end

    if all(isnan(thresholds))
        tag = sprintf('dovf_p%s', sanitizeFilenameLocal( ...
            sprintf('%.6g', percentile)));
        return;
    end

    threshold_parts = cell(1, 2);

    for k = 1:2
        if isnan(thresholds(k))
            threshold_parts{k} = sprintf('p%.6g', percentile);
        else
            threshold_parts{k} = sprintf('v%.6g', thresholds(k));
        end
    end

    tag = sanitizeFilenameLocal(sprintf( ...
        'dovf_%s_%s', threshold_parts{1}, threshold_parts{2}));
end

function [session_colors, weights, gray_colors] = ...
        buildSessionColorsLocal(group_colors, nSessions, shade_range)

    if nSessions == 1
        weights = shade_range(2);
    else
        weights = linspace(shade_range(1), shade_range(2), nSessions);
    end

    nGroups = size(group_colors, 1);
    session_colors = nan(nSessions, nGroups, 3);
    gray_colors = nan(nSessions, 3);

    for s = 1:nSessions
        w = weights(s);

        for g = 1:nGroups
            this_color = (1 - w) .* [1, 1, 1] + w .* group_colors(g, :);
            session_colors(s, g, :) = reshape(this_color, 1, 1, 3);
        end

        gray_target = [0.20, 0.20, 0.20];
        gray_colors(s, :) = (1 - w) .* [1, 1, 1] + w .* gray_target;
    end
end

function color = sessionColorLocal(session_colors, session_idx, group_idx)

    color = reshape(session_colors(session_idx, group_idx, :), 1, 3);
end

function stats = computeBasicStatsLocal(values)

    values = double(values(:));
    values = values(isfinite(values));

    stats = struct();
    stats.n = numel(values);
    stats.mean = NaN;
    stats.median = NaN;
    stats.std = NaN;
    stats.sem = NaN;
    stats.min = NaN;
    stats.max = NaN;

    if isempty(values)
        return;
    end

    stats.mean = mean(values);
    stats.median = median(values);
    stats.std = std(values, 0);
    stats.sem = stats.std ./ sqrt(numel(values));
    stats.min = min(values);
    stats.max = max(values);
end

function reg = computeLinearRegressionLocal(x, y)

    x = double(x(:));
    y = double(y(:));
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    reg = struct();
    reg.model = 'y = slope * x + intercept';
    reg.n_valid = numel(x);
    reg.is_valid = false;
    reg.reason_invalid = '';
    reg.slope = NaN;
    reg.intercept = NaN;
    reg.R2 = NaN;
    reg.SSE = NaN;
    reg.SST = NaN;
    reg.y_mean = NaN;
    reg.x_min = NaN;
    reg.x_max = NaN;
    reg.y_min = NaN;
    reg.y_max = NaN;

    if numel(x) < 2
        reg.reason_invalid = 'fewer than 2 valid points';
        return;
    end

    if numel(unique(x)) < 2
        reg.reason_invalid = 'x has fewer than 2 unique values';
        return;
    end

    p = polyfit(x, y, 1);
    yhat = p(1) .* x + p(2);
    SSE = sum((y - yhat).^2);
    SST = sum((y - mean(y)).^2);

    if SST > 0
        R2 = 1 - SSE ./ SST;
    else
        R2 = NaN;
    end

    reg.is_valid = true;
    reg.slope = p(1);
    reg.intercept = p(2);
    reg.R2 = R2;
    reg.SSE = SSE;
    reg.SST = SST;
    reg.y_mean = mean(y);
    reg.x_min = min(x);
    reg.x_max = max(x);
    reg.y_min = min(y);
    reg.y_max = max(y);
end

function central_mask = computeCentralMaskLocal(x, y, axis_info)

    x = x(:);
    y = y(:);

    if ~axis_info.use_broken_axis
        central_mask = true(size(x));
        return;
    end

    central_mask = abs(x) <= axis_info.break_start & ...
        abs(y) <= axis_info.break_start;
end

function axis_info = computeSymmetricBrokenAxisInfoLocal(x, y, options)

    use_broken_axis_display = options.use_broken_axis_display;
    break_start_prctile = options.break_start_prctile;
    broken_axis_trigger_ratio = options.broken_axis_trigger_ratio;
    tail_display_frac = options.tail_display_frac;
    break_gap_frac = options.break_gap_frac;
    min_tail_display_slope = options.min_tail_display_slope;
    force_symmetric_broken_axis = options.force_symmetric_broken_axis;

    if ~(isscalar(min_tail_display_slope) && ...
            isfinite(min_tail_display_slope) && ...
            min_tail_display_slope > 0 && min_tail_display_slope <= 1)
        error('min_tail_display_slope must be a finite scalar in (0, 1].');
    end

    if ~(isscalar(break_start_prctile) && ...
            isfinite(break_start_prctile) && ...
            break_start_prctile > 0 && break_start_prctile < 100)
        error('break_start_prctile must be a scalar between 0 and 100.');
    end

    if ~(isscalar(broken_axis_trigger_ratio) && ...
            isfinite(broken_axis_trigger_ratio) && ...
            broken_axis_trigger_ratio >= 1)
        error('broken_axis_trigger_ratio must be a finite scalar >= 1.');
    end

    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));

    axis_info = struct();
    axis_info.use_broken_axis = false;
    axis_info.break_start_prctile = break_start_prctile;
    axis_info.broken_axis_trigger_ratio = broken_axis_trigger_ratio;
    axis_info.tail_display_frac = tail_display_frac;
    axis_info.break_gap_frac = break_gap_frac;
    axis_info.min_tail_display_slope = min_tail_display_slope;
    axis_info.raw_axis_min = NaN;
    axis_info.raw_axis_max = NaN;
    axis_info.raw_abs_max = NaN;
    axis_info.break_start = NaN;
    axis_info.break_end = NaN;
    axis_info.tail_display_len = NaN;
    axis_info.break_gap = NaN;
    axis_info.plot_axis_min = NaN;
    axis_info.plot_axis_max = NaN;
    axis_info.plot_ticks = [];
    axis_info.tick_labels = {};

    if isempty(vals)
        axis_info.raw_axis_min = -1;
        axis_info.raw_axis_max = 1;
        axis_info.raw_abs_max = 1;
        axis_info.break_start = 1;
        axis_info.break_end = 1;
        axis_info.tail_display_len = 0;
        axis_info.break_gap = 0;
        axis_info.plot_axis_min = -1;
        axis_info.plot_axis_max = 1;
        [axis_info.plot_ticks, axis_info.tick_labels] = ...
            buildBrokenAxisTicksLocal(axis_info);
        return;
    end

    raw_min = min(vals);
    raw_max = max(vals);
    raw_abs_max = max(abs([raw_min, raw_max]));

    if raw_abs_max <= 0 || ~isfinite(raw_abs_max)
        raw_abs_max = 1;
    end

    if force_symmetric_broken_axis
        raw_axis_min = -raw_abs_max;
        raw_axis_max = raw_abs_max;
    else
        raw_axis_min = raw_min;
        raw_axis_max = raw_max;
    end

    abs_vals = abs(vals);
    break_start = percentileLocal(abs_vals, break_start_prctile);

    if isempty(break_start) || ~isfinite(break_start) || break_start <= 0
        break_start = raw_abs_max;
    end

    use_broken_axis = logical(use_broken_axis_display) && ...
        raw_abs_max > break_start .* broken_axis_trigger_ratio;

    axis_info.raw_axis_min = raw_axis_min;
    axis_info.raw_axis_max = raw_axis_max;
    axis_info.raw_abs_max = raw_abs_max;

    if ~use_broken_axis
        axis_info.break_start = raw_abs_max;
        axis_info.break_end = raw_abs_max;
        axis_info.tail_display_len = 0;
        axis_info.break_gap = 0;

        if raw_axis_max <= raw_axis_min
            raw_axis_min = raw_axis_min - 0.5;
            raw_axis_max = raw_axis_max + 0.5;
            axis_info.raw_axis_min = raw_axis_min;
            axis_info.raw_axis_max = raw_axis_max;
        end

        pad = 0.06 .* max(eps, raw_axis_max - raw_axis_min);
        axis_info.plot_axis_min = raw_axis_min - pad;
        axis_info.plot_axis_max = raw_axis_max + pad;
        [axis_info.plot_ticks, axis_info.tick_labels] = ...
            buildBrokenAxisTicksLocal(axis_info);
        return;
    end

    high_abs_vals = abs_vals(abs_vals > break_start);
    break_end = min(high_abs_vals);

    if isempty(break_end) || ~isfinite(break_end) || ...
            break_end <= break_start
        break_end = break_start;
    end

    raw_tail_span = max(eps, raw_abs_max - break_end);
    tail_display_len_from_frac = max(eps, break_start .* tail_display_frac);
    tail_display_len_from_slope = raw_tail_span .* min_tail_display_slope;
    tail_display_len = max( ...
        tail_display_len_from_frac, tail_display_len_from_slope);
    break_gap = max(eps, break_start .* break_gap_frac);

    axis_info.use_broken_axis = true;
    axis_info.break_start = break_start;
    axis_info.break_end = break_end;
    axis_info.tail_display_len = tail_display_len;
    axis_info.break_gap = break_gap;

    plot_min = brokenAxisTransformSignedLocal(raw_axis_min, axis_info);
    plot_max = brokenAxisTransformSignedLocal(raw_axis_max, axis_info);

    if force_symmetric_broken_axis
        plot_abs_max = max(abs([plot_min, plot_max]));
        plot_min = -plot_abs_max;
        plot_max = plot_abs_max;
    end

    pad = 0.06 .* max(eps, plot_max - plot_min);
    axis_info.plot_axis_min = plot_min - pad;
    axis_info.plot_axis_max = plot_max + pad;
    [axis_info.plot_ticks, axis_info.tick_labels] = ...
        buildBrokenAxisTicksLocal(axis_info);
end

function x_plot = brokenAxisTransformSignedLocal(x_raw, axis_info)

    x_plot = x_raw;

    if ~axis_info.use_broken_axis
        return;
    end

    sgn = sign(x_raw);
    abs_x = abs(x_raw);
    break_start = axis_info.break_start;
    break_end = axis_info.break_end;
    raw_abs_max = axis_info.raw_abs_max;
    tail_display_len = axis_info.tail_display_len;
    break_gap = axis_info.break_gap;
    x_plot_abs = abs_x;

    gap_mask = abs_x > break_start & abs_x < break_end;
    tail_mask = abs_x >= break_end;

    if break_end > break_start
        x_plot_abs(gap_mask) = break_start + ...
            (abs_x(gap_mask) - break_start) ./ ...
            (break_end - break_start) .* break_gap;
    else
        x_plot_abs(gap_mask) = break_start + break_gap;
    end

    if raw_abs_max > break_end
        x_plot_abs(tail_mask) = break_start + break_gap + ...
            (abs_x(tail_mask) - break_end) ./ ...
            (raw_abs_max - break_end) .* tail_display_len;
    else
        x_plot_abs(tail_mask) = ...
            break_start + break_gap + tail_display_len;
    end

    x_plot = sgn .* x_plot_abs;
end

function [plot_ticks, tick_labels] = buildBrokenAxisTicksLocal(axis_info)

    if axis_info.use_broken_axis
        raw_ticks = choosePrebreakTicksSignedLocal(axis_info.break_start);
        raw_ticks = uniqueWithToleranceLocal(raw_ticks, 1e-10);
        plot_ticks = brokenAxisTransformSignedLocal(raw_ticks, axis_info);
        tick_labels = composeNumericTickLabelsLocal(raw_ticks);
    else
        raw_ticks = chooseLinearTicksLocal( ...
            axis_info.raw_axis_min, axis_info.raw_axis_max);
        plot_ticks = raw_ticks;
        tick_labels = composeNumericTickLabelsLocal(raw_ticks);
    end
end

function ticks = choosePrebreakTicksSignedLocal(prebreak_abs_max)

    if isempty(prebreak_abs_max) || ~isfinite(prebreak_abs_max) || ...
            prebreak_abs_max <= 0
        ticks = 0;
        return;
    end

    rough_step = (2 .* prebreak_abs_max) ./ 4;
    step = chooseNiceNumericStepLocal(rough_step);
    max_tick = floor(prebreak_abs_max ./ step) .* step;

    if max_tick <= 0
        ticks = [-prebreak_abs_max, 0, prebreak_abs_max];
    else
        ticks = -max_tick:step:max_tick;

        if ~any(abs(ticks) < 1e-12)
            ticks = sort([ticks, 0]);
        end
    end

    ticks(abs(ticks) < 1e-12) = 0;
end

function step = chooseNiceNumericStepLocal(rough_step)

    if rough_step <= 0 || ~isfinite(rough_step)
        step = 1;
        return;
    end

    exponent = floor(log10(rough_step));
    base = 10 .^ exponent;
    candidates = [1, 2, 2.5, 5, 10] .* base;
    candidates = candidates(candidates >= rough_step);

    if isempty(candidates)
        step = 10 .* base;
    else
        step = candidates(1);
    end
end

function ticks = chooseLinearTicksLocal(axis_min, axis_max)

    if ~isfinite(axis_min) || ~isfinite(axis_max) || axis_max <= axis_min
        ticks = [-1, -0.5, 0, 0.5, 1];
        return;
    end

    rough_step = (axis_max - axis_min) ./ 4;
    step = chooseNiceNumericStepLocal(rough_step);
    tick_min = ceil(axis_min ./ step) .* step;
    tick_max = floor(axis_max ./ step) .* step;
    ticks = tick_min:step:tick_max;

    if numel(ticks) < 2
        ticks = linspace(axis_min, axis_max, 5);
    end
end

function values_unique = uniqueWithToleranceLocal(values, tolerance)

    values = values(:)';
    values_unique = [];

    for k = 1:numel(values)
        if isempty(values_unique) || ...
                all(abs(values(k) - values_unique) > tolerance)
            values_unique(end + 1) = values(k); %#ok<AGROW>
        end
    end
end

function labels = composeNumericTickLabelsLocal(values)

    labels = cell(size(values));

    for k = 1:numel(values)
        value = values(k);

        if abs(value) < 1e-12
            value = 0;
        end

        if abs(value) >= 1000 || (abs(value) > 0 && abs(value) < 0.001)
            labels{k} = sprintf('%.2g', value);
        elseif abs(value - round(value)) < 1e-10
            labels{k} = sprintf('%d', round(value));
        else
            labels{k} = sprintf('%.3g', value);
        end
    end
end

function q = percentileLocal(values, percentile)

    values = values(:);
    values = values(isfinite(values));

    if isempty(values)
        q = NaN;
        return;
    end

    values = sort(values);
    n = numel(values);

    if n == 1
        q = values;
        return;
    end

    percentile = max(0, min(100, percentile));
    position = 1 + (n - 1) .* percentile ./ 100;
    lo = floor(position);
    hi = ceil(position);

    if lo == hi
        q = values(lo);
    else
        weight = position - lo;
        q = (1 - weight) .* values(lo) + weight .* values(hi);
    end
end

function fig = plotScatterSummaryLocal(panel, allow_fit, config)

    fig = figure( ...
        'Visible', config.figure_visible, ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1, 1, config.scatter_figure_size_inches], ...
        'Renderer', 'painters', ...
        'Name', panel.key);

    ax = axes(fig);
    hold(ax, 'on');
    axis_info = panel.display_axis_info;

    drawDiagonalReferenceLocal( ...
        ax, axis_info, [0.42, 0.42, 0.42], '--', ...
        config.reference_line_width);

    zero_plot = brokenAxisTransformSignedLocal(0, axis_info);
    plot(ax, [zero_plot, zero_plot], ...
        [axis_info.plot_axis_min, axis_info.plot_axis_max], ...
        ':', 'Color', [0.72, 0.72, 0.72], ...
        'LineWidth', config.zero_line_width, ...
        'HandleVisibility', 'off');
    plot(ax, [axis_info.plot_axis_min, axis_info.plot_axis_max], ...
        [zero_plot, zero_plot], ...
        ':', 'Color', [0.72, 0.72, 0.72], ...
        'LineWidth', config.zero_line_width, ...
        'HandleVisibility', 'off');

    nSessions = size(panel.session_group, 1);
    nGroups = size(panel.session_group, 2);

    if allow_fit
        for s = 1:nSessions
            for g = 1:nGroups
                color = sessionColorLocal(config.session_colors, s, g);

                if strcmp(config.regression_mode, 'central')
                    reg = panel.session_group(s, g).regression_central;
                else
                    reg = panel.session_group(s, g).regression_all;
                end

                plotRegressionLineLocal( ...
                    ax, reg, axis_info, color, config.fit_line_width);
            end
        end
    end

    for s = 1:nSessions
        for g = 1:nGroups
            x = panel.session_group(s, g).x;
            y = panel.session_group(s, g).y;

            if isempty(x)
                continue;
            end

            x_plot = brokenAxisTransformSignedLocal(x, axis_info);
            y_plot = brokenAxisTransformSignedLocal(y, axis_info);
            color = sessionColorLocal(config.session_colors, s, g);

            drawPointsLocal( ...
                ax, x_plot, y_plot, config.scatter_point_size, ...
                color, color, config.scatter_point_alpha, 'o');
        end
    end

    axis(ax, 'square');
    xlim(ax, [axis_info.plot_axis_min, axis_info.plot_axis_max]);
    ylim(ax, [axis_info.plot_axis_min, axis_info.plot_axis_max]);
    set(ax, ...
        'XTick', axis_info.plot_ticks, ...
        'YTick', axis_info.plot_ticks, ...
        'XTickLabel', axis_info.tick_labels, ...
        'YTickLabel', axis_info.tick_labels);

    xlabel(ax, sprintf('%s (%s)', ...
        config.effect_modulation_label, panel.base_axis_label), ...
        'Interpreter', 'none', ...
        'FontName', config.font_name);
    ylabel(ax, sprintf('%s (%s)', ...
        config.effect_modulation_label, panel.comparison_axis_label), ...
        'Interpreter', 'none', ...
        'FontName', config.font_name);

    cleanAxisLocal(ax, config);

    if axis_info.use_broken_axis && config.draw_broken_axis_marks
        drawAxisBreakMarksLocal(ax, axis_info, 1.0);
    end

    if config.scatter_show_legend
        addGroupSessionLegendLocal(ax, config, 'northeastoutside');
    end

    % No title: each figure is intended to be assembled as a paper panel.
    set(findall(fig, '-property', 'FontName'), ...
        'FontName', config.font_name);
end

function drawDiagonalReferenceLocal( ...
        ax, axis_info, color, line_style, line_width)

    if ~axis_info.use_broken_axis
        raw_reference = linspace( ...
            axis_info.raw_axis_min, axis_info.raw_axis_max, 300);
        plot_reference = brokenAxisTransformSignedLocal( ...
            raw_reference, axis_info);
        plot(ax, plot_reference, plot_reference, line_style, ...
            'Color', color, ...
            'LineWidth', line_width, ...
            'HandleVisibility', 'off');
        return;
    end

    intervals = [ ...
        axis_info.raw_axis_min, min(-axis_info.break_end, axis_info.raw_axis_max); ...
        max(axis_info.raw_axis_min, -axis_info.break_start), ...
            min(axis_info.raw_axis_max, axis_info.break_start); ...
        max(axis_info.break_end, axis_info.raw_axis_min), ...
            axis_info.raw_axis_max];

    for k = 1:size(intervals, 1)
        lo = intervals(k, 1);
        hi = intervals(k, 2);

        if isfinite(lo) && isfinite(hi) && hi > lo
            raw_reference = linspace(lo, hi, 100);
            plot_reference = brokenAxisTransformSignedLocal( ...
                raw_reference, axis_info);
            plot(ax, plot_reference, plot_reference, line_style, ...
                'Color', color, ...
                'LineWidth', line_width, ...
                'HandleVisibility', 'off');
        end
    end
end

function plotRegressionLineLocal(ax, reg, axis_info, color, line_width)

    if ~reg.is_valid || ~isfinite(reg.slope) || ...
            ~isfinite(reg.intercept)
        return;
    end

    if axis_info.use_broken_axis
        raw_min = axis_info.plot_axis_min;
        raw_max = axis_info.plot_axis_max;

        if axis_info.raw_axis_min < -axis_info.break_start
            raw_min = -axis_info.break_start;
        end

        if axis_info.raw_axis_max > axis_info.break_start
            raw_max = axis_info.break_start;
        end
    else
        raw_min = axis_info.plot_axis_min;
        raw_max = axis_info.plot_axis_max;
    end

    [x_endpoints, y_endpoints] = clipLineToSquareLocal( ...
        reg.slope, reg.intercept, raw_min, raw_max);

    if numel(x_endpoints) < 2
        return;
    end

    x_plot = brokenAxisTransformSignedLocal(x_endpoints, axis_info);
    y_plot = brokenAxisTransformSignedLocal(y_endpoints, axis_info);

    plot(ax, x_plot, y_plot, '-', ...
        'Color', color, ...
        'LineWidth', line_width, ...
        'HandleVisibility', 'off');
end

function [x_endpoints, y_endpoints] = clipLineToSquareLocal( ...
        slope, intercept, raw_min, raw_max)

    x_endpoints = [];
    y_endpoints = [];

    if ~(isfinite(raw_min) && isfinite(raw_max) && raw_max > raw_min)
        return;
    end

    x_lo = raw_min;
    x_hi = raw_max;

    if abs(slope) < eps(max(1, abs(intercept)))
        if intercept < raw_min || intercept > raw_max
            return;
        end
    else
        x_at_y_min = (raw_min - intercept) ./ slope;
        x_at_y_max = (raw_max - intercept) ./ slope;
        x_from_y_lo = min(x_at_y_min, x_at_y_max);
        x_from_y_hi = max(x_at_y_min, x_at_y_max);
        x_lo = max(x_lo, x_from_y_lo);
        x_hi = min(x_hi, x_from_y_hi);
    end

    if ~(isfinite(x_lo) && isfinite(x_hi) && x_hi > x_lo)
        return;
    end

    x_endpoints = [x_lo, x_hi];
    y_endpoints = slope .* x_endpoints + intercept;

    y_endpoints = min(raw_max, max(raw_min, y_endpoints));
end

function drawAxisBreakMarksLocal(ax, axis_info, line_width)

    if ~axis_info.use_broken_axis
        return;
    end

    x_limits = ax.XLim;
    y_limits = ax.YLim;
    display_span = max(eps, max([diff(x_limits), diff(y_limits)]));
    slash_dx = display_span .* 0.010;
    slash_dy = display_span .* 0.021;
    offset = display_span .* 0.015;
    x0 = x_limits(1);
    y0 = y_limits(1);
    break_center_abs = ...
        axis_info.break_start + axis_info.break_gap ./ 2;

    x_centers = [];

    if axis_info.raw_axis_max > axis_info.break_start
        x_centers(end + 1) = break_center_abs; %#ok<AGROW>
    end

    if axis_info.raw_axis_min < -axis_info.break_start
        x_centers(end + 1) = -break_center_abs; %#ok<AGROW>
    end

    for k = 1:numel(x_centers)
        xc = x_centers(k);
        plot(ax, [xc - slash_dx, xc + slash_dx], ...
            [y0 + slash_dy, y0 - slash_dy], 'k-', ...
            'LineWidth', line_width, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off');
        plot(ax, [xc - slash_dx + offset, xc + slash_dx + offset], ...
            [y0 + slash_dy, y0 - slash_dy], 'k-', ...
            'LineWidth', line_width, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off');
    end

    y_centers = [];

    if axis_info.raw_axis_max > axis_info.break_start
        y_centers(end + 1) = break_center_abs; %#ok<AGROW>
    end

    if axis_info.raw_axis_min < -axis_info.break_start
        y_centers(end + 1) = -break_center_abs; %#ok<AGROW>
    end

    for k = 1:numel(y_centers)
        yc = y_centers(k);
        plot(ax, [x0 + slash_dx, x0 - slash_dx], ...
            [yc - slash_dy, yc + slash_dy], 'k-', ...
            'LineWidth', line_width, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off');
        plot(ax, [x0 + slash_dx, x0 - slash_dx], ...
            [yc - slash_dy + offset, yc + slash_dy + offset], 'k-', ...
            'LineWidth', line_width, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off');
    end
end

function addGroupSessionLegendLocal(ax, config, location)

    nSessions = numel(config.session_display_names);
    nGroups = numel(config.group_display_names);
    handles = gobjects(nGroups + nSessions, 1);
    labels = cell(nGroups + nSessions, 1);

    for g = 1:nGroups
        handles(g) = plot(ax, nan, nan, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', config.group_colors(g, :), ...
            'MarkerEdgeColor', config.group_colors(g, :), ...
            'LineStyle', 'none');
        labels{g} = config.group_display_names{g};
    end

    for s = 1:nSessions
        handles(nGroups + s) = plot(ax, nan, nan, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', config.session_gray_colors(s, :), ...
            'MarkerEdgeColor', config.session_gray_colors(s, :), ...
            'LineStyle', 'none');
        labels{nGroups + s} = config.session_display_names{s};
    end

    legend(ax, handles, labels, ...
        'Location', location, ...
        'Box', 'off', ...
        'Interpreter', 'none', ...
        'FontName', config.font_name, ...
        'FontSize', config.font_size);
end

function drawPointsLocal( ...
        ax, x, y, point_size, face_color, edge_color, point_alpha, marker)

    try
        scatter(ax, x, y, point_size, marker, ...
            'MarkerFaceColor', face_color, ...
            'MarkerEdgeColor', edge_color, ...
            'MarkerFaceAlpha', point_alpha, ...
            'MarkerEdgeAlpha', point_alpha, ...
            'LineWidth', 0.45, ...
            'HandleVisibility', 'off');
    catch
        scatter(ax, x, y, point_size, marker, ...
            'MarkerFaceColor', face_color, ...
            'MarkerEdgeColor', edge_color, ...
            'LineWidth', 0.45, ...
            'HandleVisibility', 'off');
    end
end

function cleanAxisLocal(ax, config)

    grid(ax, 'off');
    box(ax, 'off');
    set(ax, ...
        'TickDir', 'out', ...
        'TickLength', [0.018, 0.018], ...
        'LineWidth', config.axis_line_width, ...
        'FontName', config.font_name, ...
        'FontSize', config.font_size, ...
        'XColor', 'k', ...
        'YColor', 'k', ...
        'Layer', 'top');
end

function fig = plotSlopeSummaryLocal(slope_values, config)

    fig = figure( ...
        'Visible', config.figure_visible, ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1, 1, config.figure_size_inches], ...
        'Renderer', 'painters', ...
        'Name', 'effect slope summary');

    ax = axes(fig);
    hold(ax, 'on');

    nSessions = size(slope_values, 1);
    nComparisons = size(slope_values, 2);
    nGroups = size(slope_values, 3);
    group_offsets = centeredGroupOffsetsLocal(nGroups, config.group_offset);
    x_padding = max(0.65, max(abs(group_offsets)) + 0.35);

    y_limits = config.y_limits;

    if isempty(y_limits)
        all_values = slope_values(isfinite(slope_values));
        y_limits = automaticLimitsLocal([all_values(:); 1], true);
    end

    for c = 1:nComparisons
        for g = 1:nGroups
            x_center = config.x_positions(c) + group_offsets(g);
            values = slope_values(:, c, g);
            values = values(isfinite(values));
            jitter_by_session = deterministicJitterLocal( ...
                nSessions, config.point_jitter_width, ...
                7000 + 1000 .* c + 100 .* g);

            if isempty(values)
                continue;
            end

            switch config.plot_style
                case 'bar_points'
                    if config.show_mean
                        bar(ax, x_center, mean(values), config.bar_width, ...
                            'FaceColor', config.group_colors(g, :), ...
                            'FaceAlpha', config.bar_alpha, ...
                            'EdgeColor', 'none', ...
                            'HandleVisibility', 'off');
                    end

                case 'violin'
                    [density, y_grid] = estimateDensityLocal(values, y_limits);

                    if ~isempty(density)
                        density = density ./ max(density) .* ...
                            config.violin_width;
                        patch(ax, ...
                            [x_center - density, fliplr(x_center + density)], ...
                            [y_grid, fliplr(y_grid)], ...
                            config.group_colors(g, :), ...
                            'FaceAlpha', config.violin_alpha, ...
                            'EdgeColor', config.group_colors(g, :), ...
                            'LineWidth', 0.9, ...
                            'HandleVisibility', 'off');
                    end

                case 'points_median'
                    % No density object.
            end

            for s = 1:nSessions
                value = slope_values(s, c, g);

                if ~isfinite(value)
                    continue;
                end

                color = sessionColorLocal(config.session_colors, s, g);

                if strcmp(config.plot_style, 'bar_points')
                    edge_color = 'k';
                else
                    edge_color = 'none';
                end

                drawPointsLocal(ax, x_center + jitter_by_session(s), value, ...
                    config.point_size, color, edge_color, ...
                    config.point_alpha, 'o');
            end

            if config.show_mean && ...
                    ~strcmp(config.plot_style, 'bar_points')
                plot(ax, x_center, mean(values), 'kd', ...
                    'MarkerSize', config.mean_marker_size, ...
                    'MarkerFaceColor', 'k', ...
                    'HandleVisibility', 'off');
            end

            if config.show_median
                median_value = median(values);

                if strcmp(config.plot_style, 'bar_points')
                    plot(ax, x_center, median_value, 'ks', ...
                        'MarkerSize', 6.5, ...
                        'MarkerFaceColor', 'k', ...
                        'HandleVisibility', 'off');
                else
                    plot(ax, ...
                        [x_center - config.median_line_half_width, ...
                         x_center + config.median_line_half_width], ...
                        [median_value, median_value], '-', ...
                        'Color', config.group_colors(g, :), ...
                        'LineWidth', 2, ...
                        'HandleVisibility', 'off');
                end
            end
        end
    end

    plot(ax, ...
        [config.x_positions(1) - x_padding, ...
         config.x_positions(end) + x_padding], ...
        [1, 1], '--', ...
        'Color', [0.55, 0.55, 0.55], ...
        'LineWidth', 0.9, ...
        'HandleVisibility', 'off');

    xlim(ax, [config.x_positions(1) - x_padding, ...
              config.x_positions(end) + x_padding]);
    ylim(ax, y_limits);
    set(ax, ...
        'XTick', config.x_positions, ...
        'XTickLabel', config.comparison_labels, ...
        'TickLabelInterpreter', 'none');
    ylabel(ax, 'Regression slope', ...
        'FontName', config.font_name);
    xtickangle(ax, 20);
    cleanAxisLocal(ax, config);
    addGroupSessionLegendLocal(ax, config, 'northeastoutside');

    % No title.
    set(findall(fig, '-property', 'FontName'), ...
        'FontName', config.font_name);
end

function fig = plotDeltaSummaryLocal( ...
        delta_by_session, delta_pooled, delta_stats, config)

    fig = figure( ...
        'Visible', config.figure_visible, ...
        'Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1, 1, config.figure_size_inches], ...
        'Renderer', 'painters', ...
        'Name', 'effect delta summary');

    layout = tiledlayout(fig, 2, 1, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    nSessions = size(delta_by_session, 4);
    nComparisons = size(delta_by_session, 1);
    nGroups = size(delta_by_session, 2);
    group_offsets = centeredGroupOffsetsLocal(nGroups, config.group_offset);
    x_padding = max(0.65, max(abs(group_offsets)) + 0.35);
    axes_handles = gobjects(1, 2);

    for sign_idx = 1:2
        ax = nexttile(layout, sign_idx);
        axes_handles(sign_idx) = ax;
        hold(ax, 'on');
        overflow_info = config.overflow_info(sign_idx);
        y_limits = overflow_info.plot_y_limits;

        plot(ax, ...
            [config.x_positions(1) - x_padding, ...
             config.x_positions(end) + x_padding], ...
            [0, 0], '--', ...
            'Color', [0.60, 0.60, 0.60], ...
            'LineWidth', 0.85, ...
            'HandleVisibility', 'off');

        for c = 1:nComparisons
            for g = 1:nGroups
                x_center = config.x_positions(c) + group_offsets(g);
                values = delta_pooled{c, g, sign_idx};
                values = values(isfinite(values));
                values_plot = deltaOverflowTransformLocal( ...
                    values, overflow_info);

                if isempty(values)
                    continue;
                end

                switch config.plot_style
                    case 'bar_points'
                        if config.show_mean
                            mean_plot = deltaOverflowTransformLocal( ...
                                delta_stats(c, g, sign_idx).mean, ...
                                overflow_info);
                            bar(ax, x_center, mean_plot, ...
                                config.bar_width, ...
                                'FaceColor', config.group_colors(g, :), ...
                                'FaceAlpha', config.bar_alpha, ...
                                'EdgeColor', 'none', ...
                                'HandleVisibility', 'off');
                        end

                    case 'violin'
                        [density, y_grid] = estimateDensityLocal( ...
                            values_plot, y_limits);

                        if ~isempty(density)
                            density = density ./ max(density) .* ...
                                config.violin_width;
                            patch(ax, ...
                                [x_center - density, ...
                                 fliplr(x_center + density)], ...
                                [y_grid, fliplr(y_grid)], ...
                                config.group_colors(g, :), ...
                                'FaceAlpha', config.violin_alpha, ...
                                'EdgeColor', config.group_colors(g, :), ...
                                'LineWidth', 0.9, ...
                                'HandleVisibility', 'off');
                        end

                    case 'points_median'
                        % No density object.
                end

                for s = 1:nSessions
                    session_values = delta_by_session{c, g, sign_idx, s};
                    session_values = session_values(isfinite(session_values));
                    session_values_plot = deltaOverflowTransformLocal( ...
                        session_values, overflow_info);

                    if isempty(session_values)
                        continue;
                    end

                    jitter_seed = 10000 .* sign_idx + ...
                        1000 .* c + 100 .* g + s;
                    jitter = deterministicJitterLocal( ...
                        numel(session_values), ...
                        config.point_jitter_width, jitter_seed);
                    color = sessionColorLocal(config.session_colors, s, g);

                    if strcmp(config.plot_style, 'bar_points')
                        edge_color = 'k';
                    else
                        edge_color = 'none';
                    end

                    drawPointsLocal( ...
                        ax, x_center + jitter, session_values_plot, ...
                        config.point_size, color, edge_color, ...
                        config.point_alpha, 'o');
                end

                if config.show_mean && ...
                        ~strcmp(config.plot_style, 'bar_points')
                    mean_plot = deltaOverflowTransformLocal( ...
                        delta_stats(c, g, sign_idx).mean, overflow_info);
                    plot(ax, x_center, mean_plot, ...
                        'kd', ...
                        'MarkerSize', config.mean_marker_size, ...
                        'MarkerFaceColor', 'k', ...
                        'HandleVisibility', 'off');
                end

                if config.show_median
                    median_value = delta_stats(c, g, sign_idx).median;
                    median_plot = deltaOverflowTransformLocal( ...
                        median_value, overflow_info);

                    if strcmp(config.plot_style, 'bar_points')
                        plot(ax, x_center, median_plot, 'ks', ...
                            'MarkerSize', 6.5, ...
                            'MarkerFaceColor', 'k', ...
                            'HandleVisibility', 'off');
                    else
                        plot(ax, ...
                            [x_center - config.median_line_half_width, ...
                             x_center + config.median_line_half_width], ...
                            [median_plot, median_plot], '-', ...
                            'Color', config.group_colors(g, :), ...
                            'LineWidth', 2, ...
                            'HandleVisibility', 'off');
                    end
                end
            end
        end

        xlim(ax, [config.x_positions(1) - x_padding, ...
                  config.x_positions(end) + x_padding]);
        ylim(ax, y_limits);
        set(ax, ...
            'XTick', config.x_positions, ...
            'XTickLabel', config.comparison_labels, ...
            'TickLabelInterpreter', 'none');
        xtickangle(ax, 20);
        ylabel(ax, ['\Delta ', config.effect_modulation_label], ...
            'Interpreter', 'tex', ...
            'FontName', config.font_name);

        if strcmp(config.effect_type, 'size')
            if sign_idx == 1
                sign_label = 'Original data show suppression';
            else
                sign_label = 'Original data show facilitation';
            end
        else
            if sign_idx == 1
                sign_label = 'Original data show facilitation';
            else
                sign_label = 'Original data show suppression';
            end
        end

        sign_title = title(ax, sign_label, ...
            'Units', 'normalized', ...
            'Position', [0, 1.02, 0], ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'left', ...
            'Interpreter', 'none', ...
            'FontName', config.font_name, ...
            'FontSize', config.font_size, ...
            'FontWeight', 'normal');
        set(sign_title, 'Clipping', 'off');

        cleanAxisLocal(ax, config);

        if overflow_info.is_active
            applyDeltaOverflowTicksLocal(ax, overflow_info);
        end
    end

    legend_handles = gobjects(nGroups + nSessions, 1);
    legend_labels = cell(nGroups + nSessions, 1);

    for g = 1:nGroups
        legend_handles(g) = plot(axes_handles(1), nan, nan, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', config.group_colors(g, :), ...
            'MarkerEdgeColor', config.group_colors(g, :), ...
            'LineStyle', 'none');
        legend_labels{g} = config.group_display_names{g};
    end

    for s = 1:nSessions
        legend_handles(nGroups + s) = plot( ...
            axes_handles(1), nan, nan, 'o', ...
            'MarkerSize', 6, ...
            'MarkerFaceColor', config.session_gray_colors(s, :), ...
            'MarkerEdgeColor', config.session_gray_colors(s, :), ...
            'LineStyle', 'none');
        legend_labels{nGroups + s} = config.session_display_names{s};
    end

    lgd = legend(axes_handles(1), legend_handles, legend_labels, ...
        'Location', 'northeastoutside', ...
        'Box', 'off', ...
        'Interpreter', 'none', ...
        'FontName', config.font_name, ...
        'FontSize', config.font_size);

    try
        lgd.Layout.Tile = 'east';
    catch
        % Older MATLAB versions keep the northeastoutside position.
    end

    % Each row has a left-aligned sign-group label; there is no sgtitle.
    set(findall(fig, '-property', 'FontName'), ...
        'FontName', config.font_name);
end

function overflow_info = computeDeltaOverflowInfoLocal( ...
        delta_pooled, options)

    empty_info = struct( ...
        'is_active', false, ...
        'threshold_source', '', ...
        'threshold', NaN, ...
        'overflow_y', NaN, ...
        'plot_y_limits', [-1, 1], ...
        'raw_min', NaN, ...
        'raw_max', NaN, ...
        'n_values', 0, ...
        'n_overflow', 0);

    overflow_info = repmat(empty_info, 1, 2);

    for sign_idx = 1:2
        values = [];

        for c = 1:size(delta_pooled, 1)
            for g = 1:size(delta_pooled, 2)
                values = [values; ...
                    delta_pooled{c, g, sign_idx}(:)]; %#ok<AGROW>
            end
        end

        values = double(values(:));
        values = values(isfinite(values));
        info = empty_info;
        info.n_values = numel(values);

        if isempty(values)
            if ~isempty(options.y_limits)
                info.plot_y_limits = options.y_limits(sign_idx, :);
            end
            overflow_info(sign_idx) = info;
            continue;
        end

        info.raw_min = min(values);
        info.raw_max = max(values);

        if isfinite(options.manual_thresholds(sign_idx))
            threshold = options.manual_thresholds(sign_idx);
            info.threshold_source = 'manual';
        else
            threshold = percentileLocal(values, options.percentile);
            info.threshold_source = sprintf( ...
                'percentile_%.6g', options.percentile);
        end

        trigger_margin = (options.trigger_ratio - 1) .* ...
            max(eps, abs(threshold));
        is_active = logical(options.enabled) && isfinite(threshold) && ...
            info.raw_max > threshold + trigger_margin;

        if ~is_active
            if isempty(options.y_limits)
                info.plot_y_limits = automaticLimitsLocal(values, true);
            else
                info.plot_y_limits = options.y_limits(sign_idx, :);
            end

            info.threshold = threshold;
            overflow_info(sign_idx) = info;
            continue;
        end

        if ~isempty(options.y_limits) && ...
                options.y_limits(sign_idx, 2) < threshold
            error(['delta_y_limits row %d has an upper limit below its ' ...
                   'overflow threshold %.6g.'], sign_idx, threshold);
        end

        if ~isempty(options.y_limits) && ...
                options.y_limits(sign_idx, 1) >= threshold
            error(['delta_y_limits row %d must start below its overflow ' ...
                   'threshold %.6g.'], sign_idx, threshold);
        end

        central_values = values(values <= threshold);

        if isempty(central_values)
            central_values = threshold;
        end

        if isempty(options.y_limits)
            central_min = min([central_values(:); 0]);
            central_span = max(eps, threshold - central_min);
            lower_pad = 0.06 .* central_span;
            plot_min = central_min - lower_pad;
        else
            plot_min = options.y_limits(sign_idx, 1);
            central_span = max(eps, threshold - plot_min);
        end

        minimum_scale = max(1, max(abs([plot_min, threshold])));
        gap = max(eps .* minimum_scale, ...
            central_span .* options.gap_frac);
        top_pad = max(eps .* minimum_scale, ...
            central_span .* options.top_pad_frac);
        overflow_y = threshold + gap;

        info.is_active = true;
        info.threshold = threshold;
        info.overflow_y = overflow_y;
        info.plot_y_limits = [plot_min, overflow_y + top_pad];
        info.n_overflow = sum(values > threshold);
        overflow_info(sign_idx) = info;
    end
end

function values_plot = deltaOverflowTransformLocal(values, overflow_info)

    values_plot = values;

    if ~overflow_info.is_active
        return;
    end

    values_plot(values > overflow_info.threshold) = ...
        overflow_info.overflow_y;
end

function applyDeltaOverflowTicksLocal(ax, overflow_info)

    central_ticks = chooseLinearTicksLocal( ...
        overflow_info.plot_y_limits(1), overflow_info.threshold);
    central_ticks = central_ticks( ...
        central_ticks <= overflow_info.threshold + ...
        1e-10 .* max(1, abs(overflow_info.threshold)));

    if isempty(central_ticks)
        central_ticks = overflow_info.threshold;
    end

    ticks = [central_ticks(:)', overflow_info.overflow_y];
    labels = composeNumericTickLabelsLocal(central_ticks);
    threshold_label = formatDeltaOverflowThresholdLocal( ...
        overflow_info.threshold);
    labels{end + 1} = ['>', threshold_label];

    set(ax, 'YTick', ticks, 'YTickLabel', labels);
end

function label = formatDeltaOverflowThresholdLocal(threshold)

    % Use the threshold's own scale, not the row position, to choose the
    % displayed precision. Small thresholds retain at most two decimal
    % places; larger thresholds retain at most one. Trailing zeros and a
    % trailing decimal point are removed.
    if abs(threshold) < 1
        n_decimal_places = 2;
    else
        n_decimal_places = 1;
    end

    tolerance = 0.5 .* 10 .^ (-n_decimal_places);

    if abs(threshold) < tolerance
        threshold = 0;
    end

    label = sprintf(['%.', num2str(n_decimal_places), 'f'], threshold);
    label = regexprep(label, '0+$', '');
    label = regexprep(label, '\.$', '');

    if strcmp(label, '-0') || isempty(label)
        label = '0';
    end
end

function limits = automaticLimitsLocal(values, include_zero)

    values = double(values(:));
    values = values(isfinite(values));

    if isempty(values)
        limits = [-1, 1];
        return;
    end

    value_min = min(values);
    value_max = max(values);

    if include_zero
        value_min = min(value_min, 0);
        value_max = max(value_max, 0);
    end

    if value_max <= value_min
        pad = max(0.5, abs(value_min) .* 0.15);
    else
        pad = 0.08 .* (value_max - value_min);
    end

    limits = [value_min - pad, value_max + pad];
end

function jitter = deterministicJitterLocal(n, half_width, seed)

    if n <= 1 || half_width <= 0
        jitter = zeros(n, 1);
        return;
    end

    stream = RandStream('mt19937ar', 'Seed', seed);
    jitter = (2 .* rand(stream, n, 1) - 1) .* half_width;
end

function [density, y_grid] = estimateDensityLocal(values, y_limits)

    values = double(values(:));
    values = values(isfinite(values));
    density = [];
    y_grid = [];

    if numel(values) < 2 || max(values) == min(values)
        return;
    end

    if exist('ksdensity', 'file') == 2
        try
            [density, y_grid] = ksdensity(values, 'NumPoints', 120);
        catch
            density = [];
            y_grid = [];
        end
    end

    if isempty(density)
        nBins = min(24, max(5, round(sqrt(numel(values)))));
        [counts, edges] = histcounts(values, nBins, ...
            'Normalization', 'pdf');
        y_grid = edges(1:end - 1) + diff(edges) ./ 2;
        density = counts;
    end

    keep = isfinite(density) & isfinite(y_grid);

    if ~isempty(y_limits)
        keep = keep & y_grid >= y_limits(1) & y_grid <= y_limits(2);
    end

    density = density(keep);
    y_grid = y_grid(keep);

    if isempty(density) || max(density) <= 0
        density = [];
        y_grid = [];
        return;
    end

    density = reshape(density, 1, []);
    y_grid = reshape(y_grid, 1, []);
end

function saveFigureFilesLocal( ...
        fig, figure_files, save_fig, save_svg, save_png, png_dpi)

    drawnow;
    saved_files = {};

    if save_fig
        savefig(fig, figure_files.fig);
        saved_files{end + 1} = figure_files.fig; %#ok<AGROW>
    end

    if save_svg
        try
            exportgraphics(fig, figure_files.svg, ...
                'ContentType', 'vector', ...
                'BackgroundColor', 'none');
        catch
            print(fig, figure_files.svg, '-dsvg', '-painters');
        end
        saved_files{end + 1} = figure_files.svg; %#ok<AGROW>
    end

    if save_png
        try
            exportgraphics(fig, figure_files.png, ...
                'Resolution', png_dpi);
        catch
            print(fig, figure_files.png, '-dpng', ...
                sprintf('-r%d', png_dpi));
        end
        saved_files{end + 1} = figure_files.png; %#ok<AGROW>
    end

    fprintf('Saved figure:\n');
    for k = 1:numel(saved_files)
        fprintf('  %s\n', saved_files{k});
    end
end
