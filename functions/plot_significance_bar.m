function plot_significance_bar(time_vec, data, varargin)
% Plots a red bar above time intervals that are statistically significant.
% Also annotates the onset time of each cluster as a rounded string.
%
% Inputs:
% - time_vec: time vector corresponding to data points
% - data: p-values or boolean mask (significance)
% - cluster_size: minimum length of significant cluster (default: 1)
% - ax: axis handle (default: current axis)
% - alpha: significance threshold (default: 0.01)
% Default parameters
cluster_size = 1;
x_win = [time_vec(1), time_vec(end)];
ax = gca;
alpha = 0.01;
color = 'r';
ext_off = 0;

% Parse varargin
varStrInd = find(cellfun(@ischar, varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'cluster_size'
            cluster_size = varargin{varStrInd(iv)+1};
        case 'x_win'
            x_win = varargin{varStrInd(iv)+1};
        case 'ax'
            ax = varargin{varStrInd(iv)+1};
        case 'alpha'
            alpha = varargin{varStrInd(iv)+1};
        case 'color'
            color = varargin{varStrInd(iv)+1};

        case 'off_set'
            ext_off = varargin{varStrInd(iv)+1};
    end
end


% Threshold data if it contains p-values
is_significant = data <= alpha;

% Convert to logical transitions: 1 → significant start, -1 → end
transitions = diff([0 is_significant 0]);  % pad with zeros to catch edges

to_sig = find(transitions == 1);      % Start of significant regions
no_sig = find(transitions == -1) - 1; % End of significant regions (inclusive)

% Check for cluster size
sig_clusters = [to_sig; no_sig];  % 2 x N matrix
cluster_lengths = sig_clusters(2, :) - sig_clusters(1, :) + 1;
keep = cluster_lengths >= cluster_size;
sig_clusters = sig_clusters(:, keep);

% Get time values for start and end
sig_clusters_time = time_vec(sig_clusters);
if size(sig_clusters_time, 1) == 1
    sig_clusters_time = sig_clusters_time'; % Transform if 1 cluster present
end
% Filter clusters to be within time_vec range
time_min = x_win(1);
time_max = x_win(2);
% Find clusters that end in the time window but start before it
ends_in_window = sig_clusters_time(2,:) >= x_win(1) & sig_clusters_time(2,:) <= x_win(2);
starts_before_window = sig_clusters_time(1,:) < x_win(1);
adjust_indices = ends_in_window & starts_before_window;

% Set the start of those clusters to the start of the time window
sig_clusters_time(1, adjust_indices) = x_win(1);
keep_within_range = sig_clusters_time(1, :) >= time_min & sig_clusters_time(1, :) <= time_max; % find clusters starting in time of interest
sig_clusters_time = sig_clusters_time(:, keep_within_range);

% Plot bars on top of plot
plot_ylims = get(ax, 'YLim') + ext_off;
hold(ax, 'on');
yyaxis(ax, 'right');
set(ax, 'YTick', []);  % Remove ticks from right axis
set(ax.YAxis(2), 'Visible', 'off')
if size(sig_clusters_time, 1) == 1
    sig_clusters_time = sig_clusters_time'; % Transform if 1 cluster present
end

text_offset = 0.02 * range(plot_ylims) + ext_off; % small offset for text positioning

for i = 1:size(sig_clusters_time, 2)
    x_vals = sig_clusters_time(:, i)';
    y_vals = [plot_ylims(2), plot_ylims(2)];
    plot(ax, x_vals, y_vals, 'Color', color, 'LineWidth', 4, ...
        "LineStyle", "-", "Marker", "none");

    % Add onset time as string label (rounded to nearest whole number)
    onset_time = round(sig_clusters_time(1, i));
    label = sprintf('%d', onset_time);
    text(ax, sig_clusters_time(1, i), plot_ylims(2) + text_offset, label, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 10, 'Color', color, 'FontWeight', 'bold');
end

set(ax, 'YLim', plot_ylims); % Restore original y-limits
ax.YColor = 'k';
yyaxis(ax, 'left');
ax.YColor = 'k';
end
