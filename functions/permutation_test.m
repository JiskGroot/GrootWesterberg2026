function [p_values_cluster, sig_time_ms, sign_cluster, params] = permutation_test(data_1, data_2, params, varargin)
% permutation_test performs nonparametric cluster-based statistical testing with smoothing.
% Based on Maris & Oostenveld (2007)
%
% INPUT:
% data_1: time x trial matrix
% data_2: time x trial matrix
% params: parameter structure
% varargin: optional parameters (e.g., n_permutations, alpha_level_1, alpha_level_2, sig_time_min, time_win)
%
% OUTPUT:
% p_values_cluster: array of same length as time; NaN for non-significant times, cluster p-value for significant ones
% sig_time_ms: indices of significant time samples
% params: updated parameter structure with cluster test info

sig_time_ms = [];

% === Validate input ===
if isempty(data_1) || isempty(data_2)
    error("Ensure that you are giving two datasets as input");
end

% === Default parameters ===
n_permutations = 1000;
alpha_level_1 = 0.05; % for initial t/z thresholding
alpha_level_2 = 0.05; % for significance testing of clusters
sig_time_min = 0;     % minimum cluster length (in samples)
time_win = [1, size(data_1, 1)];
test_to_use = 'tt';   % 'tt' or 'rs'
smooth_bool = false;
nature = 0;
% === Parse varargin ===
varStrInd = find(cellfun(@ischar, varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'n_permutations'
            n_permutations = varargin{varStrInd(iv)+1};
        case 'alpha_level_1'
            alpha_level_1 = varargin{varStrInd(iv)+1};
        case 'alpha_level_2'
            alpha_level_2 = varargin{varStrInd(iv)+1};
        case 'sig_time_min'
            sig_time_min = varargin{varStrInd(iv)+1};
        case 'time_win'
            time_win = varargin{varStrInd(iv)+1};
        case 'test'
            test_to_use = varargin{varStrInd(iv)+1};
        case 'smooth'
            smooth_bool = varargin{varStrInd(iv)+1}; 
        case 'nature'
            nature = varargin{varStrInd(iv)+1};

    end
end

params.permutation_test = struct('n_permutations', n_permutations, ...
    'alpha_level_1', alpha_level_1, ...
    'alpha_level_2', alpha_level_2, ...
    'sig_time_min', sig_time_min, ...
    'time_win', time_win, ...
    'test_to_use', test_to_use);

% === Preprocessing ===
if smooth_bool
    smooth_data_1 = movmean(data_1, 15, 1, 'omitnan');
    smooth_data_2 = movmean(data_2, 15, 1, 'omitnan');
else
    smooth_data_1 = data_1;
    smooth_data_2 = data_2;
end

n_data_1_trials = size(smooth_data_1, 2);
n_data_2_trials = size(smooth_data_2, 2);
num_timepoints = size(smooth_data_1, 1);

% Check for valid number of trials
if n_data_1_trials < 5 || n_data_2_trials < 5
    error("Not enough trials in one or both conditions");
end

% Mask invalid time points (NaNs)
valid_mask = ~isnan(nanmean(smooth_data_1, 2)) & ~isnan(nanmean(smooth_data_2, 2));
smooth_data_1(~valid_mask, :) = NaN;
smooth_data_2(~valid_mask, :) = NaN;

% === Observed test statistic ===
if strcmp(test_to_use, 'tt')
    [~, ~, ~, stats] = ttest2(smooth_data_1', smooth_data_2', 'Dim', 1, ...
        'Tail', 'both', 'Vartype', 'unequal', 'Alpha', alpha_level_1);
    t_obs = stats.tstat;
elseif strcmp(test_to_use, 'rs')
    z_values = zeros(1, num_timepoints);
    for t = 1:num_timepoints
        [~, ~, stats] = ranksum(smooth_data_1(t, :), smooth_data_2(t, :), ...
            'alpha', alpha_level_1);
        z_values(t) = stats.zval;
    end
    t_obs = z_values;
else
    error('Choose a valid test: ''tt'' or ''rs''');
end

% === Cluster-forming threshold ===
t_threshold = tinv(1 - alpha_level_1 / 2, min(n_data_1_trials, n_data_2_trials) - 1);

% === Observed clusters ===
is_sig_obs = abs(t_obs) > t_threshold;
clusters_obs = bwconncomp(is_sig_obs);
n_clusters_obs = clusters_obs.NumObjects;
cluster_stats_obs = zeros(n_clusters_obs, 1);
for c = 1:n_clusters_obs
    idx = clusters_obs.PixelIdxList{c};
    cluster_stats_obs(c) = sum(abs(t_obs(idx)));
end

% === Permutation ===
combined_data = cat(2, smooth_data_1, smooth_data_2);
n_total_trials = size(combined_data, 2);
max_perm_cluster_stats = zeros(n_permutations, 1);

for perm = 1:n_permutations
    rand_idx = randperm(n_total_trials);
    perm_data_1 = combined_data(:, rand_idx(1:n_data_1_trials));
    perm_data_2 = combined_data(:, rand_idx(n_data_1_trials+1:end));

    if strcmp(test_to_use, 'tt')
        [~, ~, ~, stats_perm] = ttest2(perm_data_1', perm_data_2', 'Dim', 1, ...
            'Tail', 'both', 'Vartype', 'unequal', 'Alpha', alpha_level_2);
        t_perm = stats_perm.tstat;
    elseif strcmp(test_to_use, 'rs')
        z_values_perm = zeros(1, num_timepoints);
        for t = 1:num_timepoints
            [~, ~, stats] = ranksum(perm_data_1(t, :), perm_data_2(t, :), ...
                'alpha', alpha_level_2);
            z_values_perm(t) = stats.zval;
        end
        t_perm = z_values_perm;
    end

    is_sig_perm = abs(t_perm) > t_threshold;
    clusters_perm = bwconncomp(is_sig_perm);

    max_cluster_stat = 0;
    for c = 1:clusters_perm.NumObjects
        idx = clusters_perm.PixelIdxList{c};
        cluster_stat = sum(abs(t_perm(idx)));
        if cluster_stat > max_cluster_stat
            max_cluster_stat = cluster_stat;
        end
    end
    max_perm_cluster_stats(perm) = max_cluster_stat;
end

%% Inspect data
% Plot data
% figure
% plot(nanmean(smooth_data_1, 2)) % unprimed
% hold on
% plot(nanmean(smooth_data_2, 2))% primed
% plot(t_obs) % Plot first t-test
% yline(t_threshold)% Plot t-threshold used to identify clusters
% yline(-t_threshold)
% === Cluster-level p-values ===
cluster_pvals = ones(n_clusters_obs, 1);
for c = 1:n_clusters_obs
    cluster_pvals(c) = mean(max_perm_cluster_stats >= cluster_stats_obs(c));
end

% === Build output: p-value array for time points ===
p_values_cluster = nan(1, num_timepoints);
sign_cluster = nan(1, num_timepoints);
sig_clusters_idx = find(cluster_pvals < alpha_level_2);
for i = 1:length(sig_clusters_idx)
    idx = clusters_obs.PixelIdxList{sig_clusters_idx(i)};
    if length(idx) >= sig_time_min
        p_val = cluster_pvals(sig_clusters_idx(i));
        
        if nature == 1
            % Use sign of the mean t-stat in the cluster
            sign_dir = sign(mean(t_obs(idx), 'omitnan'));
            sign_cluster(idx) = sign_dir; 
        end
            p_values_cluster(idx) = p_val;
    end
end
sig_time_ms = find(~isnan(p_values_cluster));

% === Save cluster info to params ===
params.permutation_test.cluster_info = struct( ...
    'threshold', t_threshold);
end
