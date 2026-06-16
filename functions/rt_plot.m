function rt_plot(params)
% PLOT_RT_CDF: Plots cumulative distribution functions of reaction times
% grouped by monkey (C/H) and predictability (unpredictable: blocks 1-2,
% predictable: blocks 3-15), all individual trials pooled for plotting.
% Statistics are computed on per-session means (paired t-test + Wilcoxon).
%
% Inputs:
%   params - struct with fields:
%       .disk_dir       - path to data directory
%       .monkey_name    - 'Ca', 'He', or 'both'
%       .cut_off        - max number of blocks to include

%% Colors
% Monkey C = turquoise family, Monkey H = orange family
colors.C_unpred = [0.60 0.80 0.78];  % muted/grayish turquoise
colors.C_pred   = [0.00 0.73 0.68];  % vibrant turquoise
colors.H_unpred = [0.90 0.74 0.60];  % muted/grayish orange
colors.H_pred   = [0.95 0.45 0.05];  % vibrant orange

%% Get list of datasets
if strcmp(params.monkey_name, "both")
    file_list = dir(fullfile(params.disk_dir, '*-npy'));
elseif strcmp(params.monkey_name, "C")
    file_list = dir(fullfile(params.disk_dir, 'C*-npy'));
elseif strcmp(params.monkey_name, "H")
    file_list = dir(fullfile(params.disk_dir, 'H*-npy'));
end
dataset_names = {file_list.name};

%% Get save location
save_dir_PNG = fullfile(params.figure_save_dir, 'PNG');
save_dir_SVG = fullfile(params.figure_save_dir, 'SVG');

%% Initialize RT collectors (for plotting)
rt = struct();
rt.C_unpred = [];
rt.C_pred   = [];
rt.H_unpred = [];
rt.H_pred   = [];

%% Initialize session-level collectors (for statistics)
rt_sessions = struct('monkey', {}, 'mean_unpred', {}, 'mean_pred', {});
session_ctr = 0;

%% Main loop
for session = 1:length(dataset_names)
    dataset_name = dataset_names{session};
    dataset_path = fullfile(params.disk_dir, dataset_name);
    monkey_id    = dataset_name(1); % 'C' or 'H'

    fprintf('Processing session %d/%d: %s\n', session, length(dataset_names), dataset_name);

    %% Load required data
    try
        recording_data = matread_csv(fullfile(dataset_path, "recordinginfo.csv"));
        task_data      = matread_csv(fullfile(dataset_path, "task.csv"));
        behavior_data  = matread_npy(fullfile(dataset_path, "behavior.npy"));
    catch
        warning('Failed to load data for %s, skipping...', dataset_name);
        continue;
    end

    %% Identify valid target and distractor trials
    trials_target = ...
        (task_data.array_target_position_deg == recording_data.receptive_field_position_deg) & ...
        (behavior_data(:, 2) == 1) & ...
        ~isnan(behavior_data(:, 1));

    trials_distractor = ...
        (task_data.array_target_position_deg ~= recording_data.receptive_field_position_deg) & ...
        (task_data.catch_trial_logical == 0) & ...
        (behavior_data(:, 2) == 1) & ...
        ~isnan(behavior_data(:, 1));

    valid_trials = trials_target | trials_distractor;

    %% Split by block position (unpredictable: 1-2, predictable: 3-15)
    block_num   = task_data.block_trial_count;
    unpred_mask = valid_trials & (block_num >= 1) & (block_num <= 2);
    pred_mask   = valid_trials & (block_num >= 3) & (block_num <= min(params.cut_off, 15));

    rt_unpred = behavior_data(unpred_mask, 1);
    rt_pred   = behavior_data(pred_mask,   1);

    % Skip session if either condition has no trials
    if isempty(rt_unpred) || isempty(rt_pred)
        warning('Skipping %s: missing trials in one condition.', dataset_name);
        continue;
    end

    %% Append to plotting buckets
    if monkey_id == 'C'
        rt.C_unpred = [rt.C_unpred; rt_unpred];
        rt.C_pred   = [rt.C_pred;   rt_pred];
    elseif monkey_id == 'H'
        rt.H_unpred = [rt.H_unpred; rt_unpred];
        rt.H_pred   = [rt.H_pred;   rt_pred];
    end

    %% Collect session means for statistics
    session_ctr = session_ctr + 1;
    rt_sessions(session_ctr).monkey      = monkey_id;
    rt_sessions(session_ctr).mean_unpred = mean(rt_unpred);
    rt_sessions(session_ctr).mean_pred   = mean(rt_pred);
end

%% Compute statistics
stats = compute_rt_stats(rt_sessions);

%% Plot CDFs
fig = figure('units','normalized','outerposition',[0 0 1 1]); hold on;

plot_cdf(rt.C_unpred, colors.C_unpred, 'Monkey Ca – Unpredictable', '--');
plot_cdf(rt.C_pred,   colors.C_pred,   'Monkey Ca – Predictable',   '-');
plot_cdf(rt.H_unpred, colors.H_unpred, 'Monkey He – Unpredictable', '--');
plot_cdf(rt.H_pred,   colors.H_pred,   'Monkey He – Predictable',   '-');

xlabel('Reaction Time (ms)');
ylabel('Cumulative Probability');
title('RT Distributions by Monkey and Predictability');
legend('Location', 'southeast');
xlim([0 500]);
ylim([-.05, 1.05])
yticks([0, .5, 1])
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 13, 'Box', 'off');

%% Annotate plot with statistics
% Position annotations in top-left, stacked per monkey
annotation_x = 0.15;  % normalized figure x (adjust if needed)
y_positions  = [0.85, 0.8];  % C on top, H below
monkeys      = {'C', 'H'};
ann_colors   = {colors.C_pred, colors.H_pred};

for m = 1:2
    mk = monkeys{m};
    s  = stats.(mk);

    p_str = sprintf('W=%.0f, p=%.3f', ...
        s.W, s.p_wilcoxon);

    annotation('textbox', [annotation_x, y_positions(m), 0.5, 0.06], ...
        'String',    p_str, ...
        'Color',     ann_colors{m}, ...
        'FontSize',  11, ...
        'EdgeColor', 'none', ...
        'FontWeight','bold');
end

hold off;
saveas(fig, fullfile(save_dir_PNG, "rt_cdf.png"));
saveas(fig, fullfile(save_dir_SVG, "rt_cdf.svg"));
end

%% ========== HELPER: CDF PLOT ==========

function plot_cdf(rt_vals, color, label, linestyle)
% Plots a standard CDF (0 -> 1) for a vector of RT values

if isempty(rt_vals)
    warning('No trials found for: %s', label);
    return;
end

sorted_rt = sort(rt_vals);
cdf_y     = (1:length(sorted_rt))' / length(sorted_rt);

plot(sorted_rt, cdf_y, ...
    'Color',       color, ...
    'LineWidth',   2, ...
    'LineStyle',   linestyle, ...
    'DisplayName', label);
end

%% ========== HELPER: STATISTICS ==========

function stats = compute_rt_stats(rt_sessions)

monkeys = {'C', 'H'};
stats   = struct();

for m = 1:2
    mk   = monkeys{m};
    mask = strcmp({rt_sessions.monkey}, mk);

    unpred_means = [rt_sessions(mask).mean_unpred]';
    pred_means   = [rt_sessions(mask).mean_pred]';
    diffs        = unpred_means - pred_means;
    n            = sum(mask);

    fprintf('\nMonkey %s: n = %d sessions\n', mk, n);

    % Wilcoxon signed-rank (non-parametric)
    [p_w, ~, stats_w] = signrank(unpred_means, pred_means, 'method', 'exact');
    fprintf('  Wilcoxon signed-rank: p = %.4f, W = %.0f\n', p_w, stats_w.signedrank);

    % Effect size: matched-pairs rank-biserial correlation
    n_pairs  = sum(~isnan(diffs));
    rb       = 1 - (2 * stats_w.signedrank) / (n_pairs * (n_pairs + 1));
    fprintf('  Rank-biserial r:      %.3f\n', rb);

    % Store
    stats.(mk).n            = n;
    stats.(mk).unpred_means = unpred_means;
    stats.(mk).pred_means   = pred_means;
    stats.(mk).diffs        = diffs;
    stats.(mk).p_wilcoxon   = p_w;
    stats.(mk).W            = stats_w.signedrank;
    stats.(mk).rank_biserial = rb;
end
end