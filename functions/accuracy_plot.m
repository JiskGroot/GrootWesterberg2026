function accuracy_plot(params)
% ACCURACY_PLOT  Plots percent-correct accuracy over trials since feature change,
% separately for each monkey (C = turquoise, H = orange), with 95% CI shading.
%
% Inputs:
%   params - struct with fields:
%       .disk_dir           - path to data directory
%       .monkey_name        - 'C', 'H', or 'both'
%       .cut_off            - number of trial positions to include (e.g. 15)
%       .figure_save_dir    - root directory for saving figures

%% Colors (same family as rt_plot)
color_C       = [0.00 0.73 0.68];   % vibrant turquoise
color_C_shade = [0.60 0.80 0.78];   % muted turquoise for CI fill
color_H       = [0.95 0.45 0.05];   % vibrant orange
color_H_shade = [0.90 0.74 0.60];   % muted orange for CI fill

%% Get session list
if strcmp(params.monkey_name, "both")
    file_list = dir(fullfile(params.disk_dir, '*-npy'));
elseif strcmp(params.monkey_name, "C")
    file_list = dir(fullfile(params.disk_dir, 'C*-npy'));
elseif strcmp(params.monkey_name, "H")
    file_list = dir(fullfile(params.disk_dir, 'H*-npy'));
end
dataset_names = {file_list.name};

%% Collect per-session accuracy at each trial position
max_trials = params.cut_off;
acc_C = [];   % sessions x trial_positions (%)
acc_H = [];

for session = 1:length(dataset_names)
    dataset_name = dataset_names{session};
    dataset_path = fullfile(params.disk_dir, dataset_name);
    monkey_id    = dataset_name(1);

    fprintf('Processing session %d/%d: %s\n', session, length(dataset_names), dataset_name);

    try
        task_data     = matread_csv(fullfile(dataset_path, "task.csv"));
        behavior_data = matread_npy(fullfile(dataset_path, "behavior.npy"));
    catch
        warning('Failed to load data for %s, skipping...', dataset_name);
        continue;
    end

    % Valid trials: non-catch, non-NaN accuracy
    valid = task_data.catch_trial_logical == 0 & ~isnan(behavior_data(:, 2));

    acc_row = nan(1, max_trials);
    for t = 1:max_trials
        mask = valid & (task_data.block_trial_count == t);
        if sum(mask) > 0
            acc_row(t) = mean(behavior_data(mask, 2)) * 100;
        end
    end

    if monkey_id == 'C'
        acc_C = [acc_C; acc_row]; %#ok<AGROW>
    elseif monkey_id == 'H'
        acc_H = [acc_H; acc_row]; %#ok<AGROW>
    end
end

%% Compute session-averaged mean and 95% CI per monkey
[mean_C, lower_C, upper_C] = confidence_interval(acc_C);
[mean_H, lower_H, upper_H] = confidence_interval(acc_H);

%% Plot
x   = 1:max_trials;
fig = figure; hold on;

% Shaded confidence intervals
fill([x fliplr(x)], [upper_C fliplr(lower_C)], color_C_shade, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);
fill([x fliplr(x)], [upper_H fliplr(lower_H)], color_H_shade, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.6);

% Mean lines
line_1 = plot(x, mean_C, 'Color', color_C, 'LineWidth', 2, 'DisplayName', 'Monkey Ca');
line_2 = plot(x, mean_H, 'Color', color_H, 'LineWidth', 2, 'DisplayName', 'Monkey He');

% Chance level (1 / set size = 1/6)
line_3 = yline(100/6, '--', 'Color', [0.90 0.40 0.40], 'LineWidth', 1.5,'LineStyle', '--', 'DisplayName', 'chance');

xlabel('Trials since feature change');
ylabel('Accuracy (%)');
ylim([0 100]);
xlim([0 max_trials + 1]);
yticks([0 50 100]);
legend([line_1, line_2, line_3], {'Monkey Ca', 'Monkey He', 'chance'}, 'Location', 'southeast');
set(gca, 'TickDir', 'out', 'Box', 'off', 'FontSize', 13);
hold off;

%% Save
save_dir_PNG = fullfile(params.figure_save_dir, 'PNG');
save_dir_SVG = fullfile(params.figure_save_dir, 'SVG');
if ~exist(save_dir_PNG, 'dir'); mkdir(save_dir_PNG); end
if ~exist(save_dir_SVG, 'dir'); mkdir(save_dir_SVG); end

saveas(fig, fullfile(save_dir_PNG, 'accuracy_plot.png'));
saveas(fig, fullfile(save_dir_SVG, 'accuracy_plot.svg'));

%% Statistics — pooled across monkeys, sessions as subjects
all_acc = [acc_C; acc_H];   % n_sessions x n_trials (%)

% Repeated-measures ANOVA requires complete cases
complete = ~any(isnan(all_acc), 2);
acc_stat = all_acc(complete, :);
n_sess   = size(acc_stat, 1);

fprintf('\n=== Accuracy statistics (n = %d sessions) ===\n', n_sess);

% ---- Repeated-measures ANOVA (one within-subject factor: Trial) --------
col_names = arrayfun(@(t) sprintf('T%d', t), 1:max_trials, 'UniformOutput', false);
tbl    = array2table(acc_stat, 'VariableNames', col_names);
within = table((1:max_trials)', 'VariableNames', {'Trial'});
rm     = fitrm(tbl, sprintf('%s-%s ~ 1', col_names{1}, col_names{end}), ...
               'WithinDesign', within);
ranova_tbl = ranova(rm);

% Extract the Trial main-effect row (first row of ranova output)
F_val  = ranova_tbl.F(1);
df_num = ranova_tbl.DF(1);
df_den = ranova_tbl.DF(2);   % error row is second
p_rm   = ranova_tbl.pValue(1);

fprintf('\nRepeated-measures ANOVA — main effect of Trial:\n');
fprintf('  F(%d, %d) = %.1f,  p = %.2e\n', df_num, df_den, F_val, p_rm);

% ---- Post-hoc: paired t-test, each trial vs the next (t vs t+1) -------
fprintf('\nPost-hoc paired t-tests (consecutive trials):\n');
for t = 1:max_trials - 1
    d    = all_acc(complete, t + 1) - all_acc(complete, t);
    [~, p, ~, stats_t] = ttest(d);
    mean_diff = mean(d);
    sem_diff  = std(d) / sqrt(sum(complete));
    fprintf('  Trial %2d vs %2d:  p = %.4f,  diff = %+.2f%%  (SEM = %.2f%%),  t(%d) = %.2f\n', ...
        t, t + 1, p, mean_diff, sem_diff, stats_t.df, stats_t.tstat);
end
end
