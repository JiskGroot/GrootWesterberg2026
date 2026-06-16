% ANALYSIS.M — Complete analysis pipeline
% Reproduces all figures in Groot et al. (2026).
%
% Run main.m first to load and preprocess MUA_datamat, CSD_datamat,
% metadatamat, and params into the workspace.
%
% Figure order:
%   Figure 1C       — RT cumulative distribution functions
%   Figure 1D       — Trial since feature change accuracies
%   Figure 1E       — RT probability distribution functions (quartilized)
%   Figure 2        — Per-channel laminar MUA traces (target vs distractor)
%   Figure 3A       — Target-selection profile by RT quartile (laminar)
%   Figure 3B       — Priming effect (TE / DS) per laminar compartment
%   Figure 4B       — Neural adaptation: sequence counts per session (saved by filter_sequences_by_length)
%   Figure 4C       — Neural adaptation: MUA traces (T -> D1 -> D2 -> D3 -> D4)
%   Figure 4D       — Neural adaptation: epoch activity per session (early 58-78ms & late >78ms)
%   Figure 5A / 6A  — Laminar MUA quartile traces (target / distractor)
%   Figure 5C / 6B  — Session-average CSD maps (top) + laminar traces (bottom)

%% =========================================================================
%  SETUP  (shared parameters — run once before any figure section)
%% =========================================================================
pref_nums         = [-1, 0, 1];                  % [-1,0,1]=all | -1=non-pref | 1=pref
time_data         = params.time_data;
time_vec          = -1000:2000/2036:1000;
smooth_win        = 1;                           % sliding-mean window (samples)
sig_time_min      = 25;                          % min consecutive significant bins

if strcmp(params.rel_time_point, "ao")
    x_win = [-25 200];
elseif strcmp(params.rel_time_point, "so")
    x_win = [-75 50];
else
    error('params.rel_time_point must be "ao" or "so". Run main.m first.');
end

if any(ismember(pref_nums, 0)); pref_nums_save = 0; else; pref_nums_save = pref_nums(end); end

% Save directories
save_dir_PNG     = fullfile(params.figure_save_dir, 'MUA', 'PNG');
save_dir_SVG     = fullfile(params.figure_save_dir, 'MUA', 'SVG');
save_dir_CSD_PNG = fullfile(params.figure_save_dir, 'CSD', 'PNG');
save_dir_CSD_SVG = fullfile(params.figure_save_dir, 'CSD', 'SVG');
mkdir(fullfile(save_dir_PNG, 'lam_trace')); mkdir(fullfile(save_dir_SVG, 'lam_trace'));

% Quartile colours (target = blues, distractor = reds)
col_target = [214, 234, 248; 133, 193, 233;  52, 152, 219;  33,  97, 140] / 255;
col_dist   = [245, 200, 200; 240, 150, 130; 203,  67,  53; 120,  40,  31] / 255;

% CSD unit correction: sessions 10-29 recorded in uV, convert to nA/mm3
CSD_datamat_orig         = CSD_datamat;
idx_csd                  = metadatamat.session >= 10;
CSD_datamat(idx_csd, :)  = CSD_datamat_orig(idx_csd, :) / 1000;

% CSD delta colormap (magenta-white-green, gamma-compressed)
delta_cmap           = [1 0 1; 1 1 1; 0 1 0];
t_lin                = linspace(0, 1, 1024);
xq                   = 0.5 + 0.5 * sign(t_lin - 0.5) .* abs(2*(t_lin-0.5)).^2;
delta_cmap_interp    = interp1(linspace(0,1,3), delta_cmap, xq);

% Shared percentile bins (used by multiple figure sections)
pct_bins = [0 0.25; 0.25 0.5; 0.5 0.75; 0.75 1];

%% =========================================================================
%  FIGURE 1C — (CDFs predictable / unpredictable)
%% =========================================================================
rt_plot(params) % requires params.disk_dir, .monkey_name & .cut_off

%% =========================================================================
%  FIGURE 1D — Accuracy per trial since feature change
%% =========================================================================
accuracy_plot(params)

%% =========================================================================
%  FIGURE 1E — PDFs quartilized (predictable / unpredictable)
%% =========================================================================
prtcdf = []; urtcdf = []; urtmed = [];
quarts  = 0:0.25:1;
evalpts = 100:0.5:400;
for ii = [1:15, 17:29]
    for jj = 1:4
        prtcdf(ii,:,jj) = kde(metadatamat.rt( ...                            %#ok<SAGROW>
            quarts(jj) < metadatamat.prtprctile & ...
            metadatamat.prtprctile < quarts(jj+1) & ...
            metadatamat.session == ii), EvaluationPoints=evalpts);
        try
            urtmed(ii,jj) = nanmedian(metadatamat.rt( ...                    %#ok<SAGROW>
                quarts(jj) < metadatamat.urtprctile & ...
                metadatamat.urtprctile < quarts(jj+1) & ...
                metadatamat.session == ii));
            urtcdf(ii,:,jj) = kde(metadatamat.rt( ...                        %#ok<SAGROW>
                quarts(jj) < metadatamat.urtprctile & ...
                metadatamat.urtprctile < quarts(jj+1) & ...
                metadatamat.session == ii), EvaluationPoints=evalpts);
        catch
            urtcdf(ii,:,jj,1) = nan(length(evalpts), 1);                    %#ok<SAGROW>
        end
    end
end

fig1d = figure; nexttile; hold on;
for jj = 1:4; plot(evalpts, smoothdata(nanmean(urtcdf(:,:,jj)), 2, "movmean", 25)); end
nexttile; hold on;
for jj = 1:4; plot(evalpts, smoothdata(mean(prtcdf(:,:,jj)),   2, "movmean", 25)); end

saveas(fig1d, fullfile(save_dir_PNG, sprintf('RT_dist_upp_%s.png',      params.rel_time_point)))
saveas(fig1d, fullfile(save_dir_SVG, sprintf('RT_dist_upp_%s.svg',      params.rel_time_point)))

%% =========================================================================
%  FIGURE 2 — Per-channel laminar MUA traces
%  15-panel tiled layout per quartile x condition
%  alpha = 0.001, min sig length = 25 ms
%% =========================================================================
level_f2    = 'cat_trials';
chan_colors = [repmat([0 0 1],5,1); repmat([1 0.5 0],5,1); repmat([0 0.5 0],5,1)];
x_win_pos   = time_vec > x_win(1) & time_vec < x_win(2);
time_vals   = time_vec(x_win_pos);
for quartile = 0:4
    [primed_con, unprimed_con] = get_quartile_masks(metadatamat, quartile);
    delta_conds = {
        metadatamat.target == 1 & metadatamat.primed == 0 & unprimed_con, ...
        metadatamat.target == 0 & metadatamat.primed == 0 & unprimed_con;
        metadatamat.target == 1 & metadatamat.primed == 1 & primed_con, ...
        metadatamat.target == 0 & metadatamat.primed == 1 & primed_con;
        };

    for con = 1:2
        fig2 = figure;
        tiledlayout(15, 1, 'TileSpacing', 'none', 'Padding', 'compact');
        axs = gobjects(1, 15);
        d1 = level_of_analysis(MUA_datamat, metadatamat, level_f2, ...
            delta_conds{con,1} & ismember(metadatamat.pref, pref_nums));
        d2 = level_of_analysis(MUA_datamat, metadatamat, level_f2, ...
            delta_conds{con,2} & ismember(metadatamat.pref, pref_nums));

        for i = 1:15
            axs(i) = nexttile; hold on;
            y1 = nanmean(d1(i, x_win_pos, :), 3);
            y2 = nanmean(d2(i, x_win_pos, :), 3);
            plot(time_vals, y1, 'Color', [0 0 0],       'LineWidth', 1.2);
            plot(time_vals, y2, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
            dmask = diff([0, y1 > y2, 0]);
            startIdx = find(dmask == 1); endIdx = find(dmask == -1) - 1;
            for k = 1:length(startIdx)
                ir = startIdx(k):endIdx(k);
                fill([time_vals(ir) fliplr(time_vals(ir))], [y1(ir) fliplr(y2(ir))], ...
                    chan_colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
            end
            p = running_rs(squeeze(d1(i,:,:))', squeeze(d2(i,:,:))', 0.01);
            plot_significance_bar(time_vec, p, 'cluster_size', 25, ...
                'x_win', x_win, 'ax', axs(i), 0.01);

            set(axs(i), 'XLim', x_win, 'YLim', [-.25, 2.5], ...
                'YTick', [], 'YColor', 'none', 'Box', 'off');

            if i < 15
                axs(i).XAxis.Visible = 'off';
            else
                set(axs(i), 'TickDir', 'out', 'XTick', [0, 100, 200]);
                xlabel('Time from array onset (ms)');
            end
        end

        linkaxes(axs, 'xy');

        % Add 1 a.u. scale bar to top panel on last quartile - predictable
        if con == 2 && quartile == 4
            axes(axs(1));
            scale_x      = x_win(2) * 0.97;
            scale_y_bot  = 1.5;
            scale_height = 1;
            plot([scale_x scale_x], [scale_y_bot, scale_y_bot + scale_height], ...
                'k-', 'LineWidth', 1.5);
            text(scale_x + 4, scale_y_bot + scale_height/2, '1 a.u.', ...
                'FontSize', 7, 'VerticalAlignment', 'middle');
        end

        saveas(fig2, fullfile(save_dir_PNG, 'lam_trace', ...
            sprintf('trace_%s_%d_%d_%d.png', params.rel_time_point, quartile, con, pref_nums_save)))
        saveas(fig2, fullfile(save_dir_SVG, 'lam_trace', ...
            sprintf('trace_%s_%d_%d_%d.svg', params.rel_time_point, quartile, con, pref_nums_save)))
    end
end

%% =========================================================================
%  FIGURE 2 — TST statistics: unpredictable vs predictable per compartment
%  Per-session onset detected as first timepoint where (target - distractor)
%  exceeds baseline mean+2SD for >= sig_time_min consecutive bins.
%  Paired two-sided t-test (alpha = 0.01) for "All" trials and each quartile.
%% =========================================================================
compartments_tst = {"upper", "middle", "lower"};
alpha_tst        = 0.01;
bsl_win_tst      = time_vec >= -250 & time_vec < 0;

fprintf('\n=== TST: unpredictable vs predictable (paired signedrank, alpha=%.3f) ===\n', alpha_tst);

for quartile = 0:4
    [primed_con_q, unprimed_con_q] = get_quartile_masks(metadatamat, quartile);
    if quartile == 0
        q_label = 'All';
    else
        q_label = sprintf('Q%d', quartile);
    end
    fprintf('\n  [%s]\n', q_label);

    for c_idx = 1:numel(compartments_tst)
        comp = compartments_tst{c_idx};

        d_tgt_unp = level_of_analysis(MUA_datamat, metadatamat, comp, ...
            metadatamat.target == 1 & metadatamat.primed == 0 & unprimed_con_q & ...
            ismember(metadatamat.pref, pref_nums));   % [sessions x time]
        d_dst_unp = level_of_analysis(MUA_datamat, metadatamat, comp, ...
            metadatamat.target == 0 & metadatamat.primed == 0 & unprimed_con_q & ...
            ismember(metadatamat.pref, pref_nums));
        d_tgt_pr  = level_of_analysis(MUA_datamat, metadatamat, comp, ...
            metadatamat.target == 1 & metadatamat.primed == 1 & primed_con_q & ...
            ismember(metadatamat.pref, pref_nums));
        d_dst_pr  = level_of_analysis(MUA_datamat, metadatamat, comp, ...
            metadatamat.target == 0 & metadatamat.primed == 1 & primed_con_q & ...
            ismember(metadatamat.pref, pref_nums));

        diff_unp  = d_tgt_unp - d_dst_unp;   % [sessions x time]
        diff_pr   = d_tgt_pr  - d_dst_pr;

        nSess_tst = size(diff_unp, 1);
        tst_unp   = nan(nSess_tst, 1);
        tst_pr    = nan(nSess_tst, 1);

        for s = 1:nSess_tst
            tst_unp(s) = compute_tst_onset(diff_unp(s,:), time_vec, bsl_win_tst, x_win_pos, sig_time_min);
            tst_pr(s)  = compute_tst_onset(diff_pr(s,:),  time_vec, bsl_win_tst, x_win_pos, sig_time_min);
        end

        valid     = ~isnan(tst_unp) & ~isnan(tst_pr);
        tst_unp_v = tst_unp(valid);
        tst_pr_v  = tst_pr(valid);

        [p_val, ~, stats_t] = signrank(tst_unp_v, tst_pr_v, 'tail', 'both');

        sig_marker = '';
        if p_val < alpha_tst; sig_marker = ' *'; end
        fprintf('  %-8s  W=%6.1f  p=%7.4f%s\n', ...
            comp, stats_t.signedrank, p_val, sig_marker);
    end
end

%% =========================================================================
%  FIGURE 3A — Target-selection profile per RT quartile
%  6 panels: 3 compartments x 2 primed conditions
%  alpha = 0.01, min. clus. len. = 25 ms
%% =========================================================================
levels_3a       = ["upper", "middle", "lower"];
primed_conds_3a = [0, 1];
line_colors     = [0.8 0.8 0.8; 0.6 0.6 0.6; 0.4 0.4 0.4; 0.2 0.2 0.2];
ax_handles_3a   = gobjects(length(levels_3a) * length(primed_conds_3a), 1);
ax_idx          = 1;

fig3a = figure;
for lvl = 1:length(levels_3a)
    for p_idx = 1:length(primed_conds_3a)
        ax     = nexttile; hold on;
        level  = levels_3a(lvl);
        primed = primed_conds_3a(p_idx);
        if primed == 0; pct_field = "urtprctile"; else; pct_field = "prtprctile"; end

        data_sets = cell(1, size(pct_bins, 1));
        for i = 1:size(pct_bins, 1)
            lo = pct_bins(i,1); hi = pct_bins(i,2);
            if lo == 0; pct_mask = metadatamat.(pct_field) < hi;
            else;       pct_mask = metadatamat.(pct_field) > lo & metadatamat.(pct_field) < hi;
            end
            mt = metadatamat.target == 1 & metadatamat.primed == primed & ...
                pct_mask & ismember(metadatamat.pref, pref_nums);
            md = metadatamat.target == 0 & metadatamat.primed == primed & ...
                pct_mask & ismember(metadatamat.pref, pref_nums);
            data_sets{i} = level_of_analysis(MUA_datamat, metadatamat, level, mt) ...
                         - level_of_analysis(MUA_datamat, metadatamat, level, md);
        end
        kwp = running_kwt(data_sets);
        plot_significance_bar(time_vec, kwp, 'cluster_size', sig_time_min, 'x_win', x_win,'alpha', 0.01, 'ax', ax)
        for i = 1:length(data_sets)
            [md, ~, ~] = confidence_interval(data_sets{i});
            plot(time_vec, smoothdata(md, 'movmean', smooth_win), ...
                'Color', line_colors(i,:), 'LineWidth', 2, 'LineStyle', '-')
            xticks([0, 100, 200]);
            yticks([0, .25, .5, .75]);

        end
        yline(0);
        ax_handles_3a(ax_idx) = ax; ax_idx = ax_idx + 1;
    end
end
linkaxes(ax_handles_3a, 'xy')
set(ax_handles_3a, 'TickDir', 'out', 'Box', 'off')
set(ax_handles_3a(1), 'XLim', x_win)
set(ax_handles_3a(5), 'YLim', [-0.2 0.75])
% set(ax_handles_3a(5), 'Ytick', [0, .25, .5, .75])
% set(ax_handles_3a(1), 'Xtick', [0, 100, 200])
ax_handles_3a(5).XLabel.String = 'Time from Array Onset (ms)';
ax_handles_3a(5).YLabel.String = 'Target Selection Profile (au)';
saveas(fig3a, fullfile(save_dir_PNG, sprintf('T_D_delta_q_lam_%s_%d.png', params.rel_time_point, pref_nums_save)))
saveas(fig3a, fullfile(save_dir_SVG, sprintf('T_D_delta_q_lam_%s_%d.svg', params.rel_time_point, pref_nums_save)))


%% =========================================================================
%  FIGURE 3B — Priming effect (TE / DS) per laminar compartment
%  4 panels: upper / middle / lower / session-level
%  Plots (primed - unprimed) for both target and distractor
%% =========================================================================
fig_3b    = figure;
levels_3b = {"upper", "middle", "lower", "session"};
ax_3b     = gobjects(1, numel(levels_3b));

for i = 1:numel(levels_3b)
    level    = levels_3b{i};
    ax_3b(i) = nexttile; hold on;

    data1 = level_of_analysis(MUA_datamat, metadatamat, level, ...
        metadatamat.target == 1 & metadatamat.primed == 1 & ismember(metadatamat.pref, pref_nums)) ...
        - level_of_analysis(MUA_datamat, metadatamat, level, ...
        metadatamat.target == 1 & metadatamat.primed == 0 & ismember(metadatamat.pref, pref_nums));
    data2 = level_of_analysis(MUA_datamat, metadatamat, level, ...
        metadatamat.target == 0 & metadatamat.primed == 1 & ismember(metadatamat.pref, pref_nums)) ...
        - level_of_analysis(MUA_datamat, metadatamat, level, ...
        metadatamat.target == 0 & metadatamat.primed == 0 & ismember(metadatamat.pref, pref_nums));

    [m1, l1, u1] = confidence_interval(data1);
    [m2, l2, u2] = confidence_interval(data2);
    bsl1 = nanmean(data1(:, time_data < 0 & time_data > -250), 2);
    bsl2 = nanmean(data2(:, time_data < 0 & time_data > -250), 2);
    p1 = running_rs_null(data1, bsl1, 0.01);
    p2 = running_rs_null(data2, bsl2, 0.01);

    plot_significance_bar(time_vec, p1, 'cluster_size', 20, ...
        'x_win', x_win, 'ax', ax_3b(i), 'color', 'b', 'alpha', 0.01);
    plot_significance_bar(time_vec, p2, 'cluster_size', 20, ...
        'x_win', x_win, 'ax', ax_3b(i), 'color', 'r', 'alpha', 0.01);
    plot(time_vec, l1, 'Color','b', 'LineStyle','-'); plot(time_vec, u1, 'Color','b', 'LineStyle','-');
    plot(time_vec, l2, 'Color','r', 'LineStyle','-'); plot(time_vec, u2, 'Color','r', 'LineStyle','-');
    yline(0, '--', 'Color', [0.5 0.5 0.5]);
    plot(time_vec, smoothdata(m1, 'movmean', 1), ...
        'Color', col_target(4,:), 'LineWidth', 1.5, 'LineStyle','-', 'Marker','none');
    plot(time_vec, smoothdata(m2, 'movmean', 1), ...
        'Color', col_dist(4,:), 'LineWidth', 1.5, 'LineStyle','-', 'Marker','none');
    xlabel('Time from array onset (ms)'); ylabel('Primed - Unprimed (a.u.)');
end
linkaxes(ax_3b, 'xy')
set(ax_3b, 'XLim', x_win, 'YLim', [-.2 .2], 'TickDir', 'out')
saveas(fig_3b, fullfile(save_dir_PNG, sprintf('TE_DS_%s_%d.png', params.rel_time_point, pref_nums_save)))
saveas(fig_3b, fullfile(save_dir_SVG, sprintf('TE_DS_%s_%d.svg', params.rel_time_point, pref_nums_save)))


%% =========================================================================
%  FIGURE 4B / 4C / 4D — Neural adaptation sequence analysis
%  4B: session count histogram (saved inside filter_sequences_by_length)
%  4C: group-average MUA traces per compartment
%  4D: epoch-average activity per session (early 58-78ms and late >78ms)
%% =========================================================================

% Extract T -> D -> D -> D -> D sequences from every session
allSeqMeta = cell(1, 29);
for s = 1:29
    allSeqMeta{s} = extract_target_distractor_sequences(metadatamat, s);
end

nTrials     = 5;   % T, D1, D2, D3, D4
nOccurences = 15;  % minimum sequences required per session

% filter_sequences_by_length also creates and saves Figure 4B internally
allSeqMeta_filt = filter_sequences_by_length(allSeqMeta, nTrials, nOccurences, ...
    save_dir_PNG, save_dir_SVG);

% ---- Figure 4C: group-average traces per laminar compartment ----------
nChannels    = 15;
nSessions_ad = 29;
chan_data_seq = cell(1, nChannels);

for n_chan = 1:nChannels
    avg_data = nan(nSessions_ad, length(time_data), nTrials);
    for n_trial = 1:nTrials
        for i_session = 1:size(allSeqMeta_filt, 2)
            if ~iscell(allSeqMeta_filt{1, i_session}), continue; end
            indices   = [];
            data_sess = allSeqMeta_filt{1, i_session};
            for i_occ = 1:size(data_sess, 2)
                data_occ = data_sess{1, i_occ};
                trials   = unique(data_occ.trial);
                idx_t    = find(metadatamat.session == i_session & ...
                    metadatamat.trial == trials(n_trial+1) & ...
                    metadatamat.channel == n_chan);
                indices  = [indices, idx_t]; %#ok<AGROW>
            end
            avg_data(i_session, :, n_trial) = nanmean(MUA_datamat(indices, :), 1);
        end
    end
    chan_data_seq{n_chan} = avg_data;
end

seq_colors          = zeros(nTrials, 3);
seq_colors(1,:)     = [0 0 1];
steps               = linspace(0.3, 0.8, nTrials-1)';
seq_colors(2:end,:) = repmat(steps, 1, 3);

for compartment = ["upper", "middle", "lower"]
    switch compartment
        case 'upper'; chan_int = 1:5;
        case 'middle'; chan_int = 6:10;
        case 'lower'; chan_int = 11:15;
    end
    nChans_ad = numel(chan_int);
    data_avg  = nan(nTrials, length(time_data));
    for n_trial = 1:nTrials
        data_temp = nan(nChans_ad * nSessions_ad, length(time_data));
        ctr = 1;
        for s = 1:nSessions_ad
            for ic = 1:nChans_ad
                data_temp(ctr,:) = squeeze(chan_data_seq{chan_int(ic)}(s, :, n_trial));
                ctr = ctr + 1;
            end
        end
        [data_avg(n_trial,:), ~, ~] = confidence_interval(data_temp);
    end
    fig = figure; hold on;
    for i_trial = 1:nTrials
        if i_trial == 1
            plot(time_data, data_avg(i_trial,:), 'Color', seq_colors(i_trial,:), 'LineWidth', 1, 'LineStyle','--');
        else
            plot(time_data, data_avg(i_trial,:), 'Color', seq_colors(i_trial,:), 'LineWidth', 1, 'LineStyle','-');
        end
    end
    xlim([-25, 200]); ylim([0, 1.5]); yticks([0, 0.5, 1, 1.5])
    xticks([-25, 0, 50, 100, 150, 200])
    xlabel('Time (ms)'); ylabel('MUA (a.u.)');
    title(sprintf('Group average (%s channels)', compartment));
    set(gca, 'TickDir', 'out'); hold off
    saveas(fig, fullfile(save_dir_PNG, sprintf('neural_adap_%d_%s.png', nTrials, compartment)));
    saveas(fig, fullfile(save_dir_SVG, sprintf('neural_adap_%d_%s.svg', nTrials, compartment)));

    if strcmp(compartment, 'middle')
        fig_in = figure; hold on;
        for i_trial = 1:nTrials
            if i_trial == 1
                plot(time_data, data_avg(i_trial,:), 'Color', seq_colors(i_trial,:), 'LineWidth', 1, 'LineStyle', '--');
            else
                plot(time_data, data_avg(i_trial,:), 'Color', seq_colors(i_trial,:), 'LineWidth', 1, 'LineStyle', '-');
            end
        end
        xlim([50, 75]); ylim([1.25, 1.5]); yticks([1.25, 1.5]); xticks([50, 75])
        xlabel('Time (ms)'); ylabel('MUA (a.u.)'); set(gca, 'TickDir', 'out'); hold off
        saveas(fig_in, fullfile(save_dir_PNG, sprintf('neural_adap_%d_%s_inset.png', nTrials, compartment)));
        saveas(fig_in, fullfile(save_dir_SVG, sprintf('neural_adap_%d_%s_inset.svg', nTrials, compartment)));
    end
end

% ---- Figure 4D: epoch activity — early (58-78ms) and late (>78ms) -----
% Compartment info used by both epoch loops
comp_info_ad = struct('name', {"upper","middle","lower"}, ...
                      'chans', {1:5, 6:10, 11:15}, 'idx', {1, 2, 3});

time_epoch_early = time_data >= 58 & time_data <= 78;   % granular onset window
time_epoch_late  = time_data >= 78 & time_data <= 200;  % sustained window

chan_data_early = cell(nSessions_ad, 3, nTrials);
chan_data_late  = cell(nSessions_ad, 3, nTrials);

for epoch_idx = 1:2
    if epoch_idx == 1
        time_epoch  = time_epoch_early;
        chan_data_ep = chan_data_early;
    else
        time_epoch  = time_epoch_late;
        chan_data_ep = chan_data_late;
    end
    for c_info = comp_info_ad
        for n_trial = 1:nTrials
            for i_session = 1:size(allSeqMeta_filt, 2)
                if ~iscell(allSeqMeta_filt{1, i_session}), continue; end
                indices   = [];
                data_sess = allSeqMeta_filt{1, i_session};
                for i_occ = 1:size(data_sess, 2)
                    data_occ = data_sess{1, i_occ};
                    trials   = unique(data_occ.trial);
                    idx_t    = find(metadatamat.session == i_session & ...
                        metadatamat.trial == trials(n_trial+1) & ...
                        ismember(metadatamat.channel, c_info.chans));
                    indices  = [indices; idx_t]; %#ok<AGROW>
                end
                chan_data_ep{i_session, c_info.idx, n_trial} = ...
                    nanmean(MUA_datamat(indices, time_epoch), 2);
            end
        end
    end
    if epoch_idx == 1
        chan_data_early = chan_data_ep;
    else
        chan_data_late  = chan_data_ep;
    end
end

% Plot late-epoch session boxplots (saves internally when save_bool=1)
plot_chan_data(chan_data_late, 'late', 'session', 0, 1, save_dir_PNG, save_dir_SVG, 0);
plot_chan_data(chan_data_early, 'early','session', 0, 1, save_dir_PNG, save_dir_SVG, 0);

% Plot early-vs-late session difference (saves internally when save_bool=1)
plot_chan_diff_session(chan_data_early, chan_data_late, 1, 1, save_dir_PNG, save_dir_SVG);


%% =========================================================================
%  FIGURE 5A / 6A — Laminar MUA traces by RT quartile
%  12 panels: 3 compartments x 4 conditions (unprimed-T, primed-T, unprimed-D, primed-D)
%% =========================================================================
fig_5a6a = figure;
ax_5a6a  = gobjects(12, 1);
ctrlam   = 0;

for ii = 1:3
    switch ii; case 1; level = "upper"; case 2; level = "middle"; case 3; level = "lower"; end
    for pair = 1:4   % 1=unprimed-T  2=primed-T  3=unprimed-D  4=primed-D
        ctrlam          = ctrlam + 1;
        ax_5a6a(ctrlam) = nexttile; hold on;
        if     pair == 1; is_target = 1; is_primed = 0; pct_f = "urtprctile"; c_set = col_target;
        elseif pair == 2; is_target = 1; is_primed = 1; pct_f = "prtprctile"; c_set = col_target;
        elseif pair == 3; is_target = 0; is_primed = 0; pct_f = "urtprctile"; c_set = col_dist;
        else;             is_target = 0; is_primed = 1; pct_f = "prtprctile"; c_set = col_dist;
        end
        data_q = cell(1, 4);
        for iq = 1:4
            lo = pct_bins(iq,1); hi = pct_bins(iq,2);
            if lo == 0; pm = metadatamat.(pct_f) < hi;
            else;       pm = metadatamat.(pct_f) > lo & metadatamat.(pct_f) < hi;
            end
            data_q{iq} = level_of_analysis(MUA_datamat, metadatamat, level, ...
                metadatamat.target == is_target & metadatamat.primed == is_primed & ...
                 pm & ismember(metadatamat.pref, pref_nums));
        end
        kwp = running_kwt(data_q);
        plot_significance_bar(time_vec, kwp, 'cluster_size', sig_time_min, ...
            'x_win', x_win, 'ax', ax_5a6a(ctrlam))
        for iq = 1:4
            [m, ~, ~] = confidence_interval(data_q{iq});
            plot(time_vec, smoothdata(m, 'movmean', smooth_win), 'Color', c_set(iq,:), 'LineWidth', 2, 'LineStyle','-')
        end
        yline(0)
    end
end
linkaxes(ax_5a6a(1:4),  'xy'); set(ax_5a6a(1), 'XLim', x_win, 'YLim', [0 1.5])
linkaxes(ax_5a6a(5:8),  'xy'); set(ax_5a6a(5), 'XLim', x_win, 'YLim', [0 2.25])
linkaxes(ax_5a6a(9:12), 'xy'); set(ax_5a6a(9), 'XLim', x_win, 'YLim', [0 1.5])
set(ax_5a6a, 'TickDir', 'out')
saveas(fig_5a6a, fullfile(save_dir_PNG, sprintf('T_D_u_p_lam_%s_%d.png', params.rel_time_point, pref_nums_save)))
saveas(fig_5a6a, fullfile(save_dir_SVG, sprintf('T_D_u_p_lam_%s_%d.svg', params.rel_time_point, pref_nums_save)))


%% =========================================================================
%  FIGURE 5C / 6B (top) — Session-average CSD maps per quartile
%% =========================================================================
sig_time_min_csd = 10;
time_pos_csd    = time_vec > x_win(1) & time_vec < x_win(2);
time_vals_csd   = time_vec(time_pos_csd);
level_csd       = 'session2d';
titles_csd      = {'Target unprimed', 'Target primed', 'Dist unprimed', 'Dist primed'};

for quartile = 0:4
    [primed_con, unprimed_con] = get_quartile_masks(metadatamat, quartile);
    nd_conds = {
        metadatamat.target == 1 & metadatamat.primed == 0 & unprimed_con;
        metadatamat.target == 1 & metadatamat.primed == 1 & primed_con;
        metadatamat.target == 0 & metadatamat.primed == 0 & unprimed_con;
        metadatamat.target == 0 & metadatamat.primed == 1  & primed_con;
        };
    all_data_csd = cell(1, 4);
    for i = 1:4
        dt = level_of_analysis(CSD_datamat, metadatamat, level_csd, ...
            nd_conds{i} & ismember(metadatamat.pref, pref_nums));
        dt = squeeze(nanmean(dt, 3));
        all_data_csd{i} = SMOOTH_2D(dt);
        all_data_csd{i} = all_data_csd{i}(:, time_pos_csd);
    end
    fig5b6c = figure;
    axs_csd = gobjects(1, 4);
    for i = 1:4
        axs_csd(i) = nexttile; hold on;
        imagesc(time_vals_csd, 1:size(all_data_csd{i},1), all_data_csd{i}, [-225, 225]);
        colormap(axs_csd(i), flipud(jet)); colorbar(axs_csd(i), 'eastoutside');
        title(titles_csd{i});
    end
    linkaxes(axs_csd, 'xy');
    for ax = axs_csd
        set(ax, 'YDir', 'reverse', 'XLim', x_win, 'TickDir', 'out', 'YTickLabel', []);
    end
    saveas(fig5b6c, fullfile(save_dir_CSD_PNG, sprintf('CSD_%d_%s_%d.png', quartile, params.rel_time_point, pref_nums_save)))
    saveas(fig5b6c, fullfile(save_dir_CSD_SVG, sprintf('CSD_%d_%s_%d.svg', quartile, params.rel_time_point, pref_nums_save)))
end


%% =========================================================================
%  FIGURE 5C / 6B (bottom) — CSD laminar time series
%  Significance against initial response window (48-78 ms)
%% =========================================================================
alpha_csd    = 0.001;
time_initial = time_vec >= 48 & time_vec < 78;
time_integ   = [78, 200];

for target = 0:1
    target_con = metadatamat.target == target;
    for quartile = 0:4
        [primed_con, unprimed_con] = get_quartile_masks(metadatamat, quartile);

        d_upp_up = level_of_analysis(CSD_datamat, metadatamat, "upper",  unprimed_con & target_con);
        d_upp_p  = level_of_analysis(CSD_datamat, metadatamat, "upper",  primed_con   & target_con);
        d_mid_up = level_of_analysis(CSD_datamat, metadatamat, "middle", unprimed_con & target_con);
        d_mid_p  = level_of_analysis(CSD_datamat, metadatamat, "middle", primed_con   & target_con);
        d_low_up = level_of_analysis(CSD_datamat, metadatamat, "lower",  unprimed_con & target_con);
        d_low_p  = level_of_analysis(CSD_datamat, metadatamat, "lower",  primed_con   & target_con);

        ini = @(d) reshape(d(:, time_initial), 1, []);
        p1_up = running_rs_null(d_upp_up, ini(d_upp_up), alpha_csd);
        p2_up = running_rs_null(d_mid_up, ini(d_mid_up), alpha_csd);
        p3_up = running_rs_null(d_low_up, ini(d_low_up), alpha_csd);
        p1_p  = running_rs_null(d_upp_p,  ini(d_upp_p),  alpha_csd);
        p2_p  = running_rs_null(d_mid_p,  ini(d_mid_p),  alpha_csd);
        p3_p  = running_rs_null(d_low_p,  ini(d_low_p),  alpha_csd);

        fig5b6c2 = figure; hold on
        plot(time_vec, nanmean(d_upp_up,1), 'Color', [198 198 229]/255)
        plot(time_vec, nanmean(d_mid_up,1), 'Color', [229 226 198]/255)
        plot(time_vec, nanmean(d_low_up,1), 'Color', [207 229 198]/255)
        plot(time_vec, nanmean(d_upp_p, 1), 'Color', [ 51   0 255]/255)
        plot(time_vec, nanmean(d_mid_p, 1), 'Color', [255 128   0]/255)
        plot(time_vec, nanmean(d_low_p, 1), 'Color', [ 60 255   0]/255)
        xlim([-25 200]); xline(0)
        ylim([-200, 100]); yticks([-200 -150 -100 -50 0 50 100]); xticks([0 100 200])
        plot_significance_bar(time_vec, p1_up, 'cluster_size', sig_time_min_csd, ...
            'color', [198 198 229]/255, 'x_win', time_integ, 'off_set', -8);
        plot_significance_bar(time_vec, p2_up, 'cluster_size', sig_time_min_csd, ...
            'color', [229 226 198]/255, 'x_win', time_integ, 'off_set', -7);
        plot_significance_bar(time_vec, p3_up, 'cluster_size', sig_time_min_csd, ...
            'color', [207 229 198]/255, 'x_win', time_integ, 'off_set', -6);
        plot_significance_bar(time_vec, p1_p,  'cluster_size', sig_time_min_csd, ...
            'color', [ 51   0 255]/255, 'x_win', time_integ, 'off_set', -4);
        plot_significance_bar(time_vec, p2_p,  'cluster_size', sig_time_min_csd, ...
            'color', [255 128   0]/255, 'x_win', time_integ, 'off_set', -3);
        plot_significance_bar(time_vec, p3_p,  'cluster_size', sig_time_min_csd, ...
            'color', [ 60 255   0]/255, 'x_win', time_integ, 'off_set', -2);
        set(gca, 'YDir', 'reverse', 'TickDir', 'out'); hold off
        saveas(fig5b6c2, fullfile(save_dir_CSD_PNG, ...
            sprintf('lat_%d_%d_CSD_%s.png', target, quartile, params.rel_time_point)))
        saveas(fig5b6c2, fullfile(save_dir_CSD_SVG, ...
            sprintf('lat_%d_%d_CSD_%s.svg', target, quartile, params.rel_time_point)))
    end
end

%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================
function tst = compute_tst_onset(diff_trace, time_vec, bsl_win, win_mask, min_dur)
% Returns the first time (ms) within win_mask where diff_trace exceeds
% baseline mean + 2 SD for at least min_dur consecutive samples. NaN if absent.
    bsl       = diff_trace(bsl_win);
    threshold = mean(bsl, 'omitnan') + 2 * std(bsl, 'omitnan');
    win_diff  = diff_trace(win_mask);
    win_time  = time_vec(win_mask);
    above     = win_diff > threshold & win_diff > 0;
    tst       = NaN;
    for t = 1:(numel(above) - min_dur + 1)
        if all(above(t : t + min_dur - 1))
            tst = win_time(t);
            return;
        end
    end
end

function [primed_con, unprimed_con] = get_quartile_masks(metadatamat, quartile)
switch quartile
    case 0
        primed_con   = true(size(metadatamat,1),1);
        unprimed_con = true(size(metadatamat,1),1);
    case 1
        primed_con   = metadatamat.prtprctile < 0.25;
        unprimed_con = metadatamat.urtprctile < 0.25;
    case 2
        primed_con   = metadatamat.prtprctile >= 0.25 & metadatamat.prtprctile < 0.5;
        unprimed_con = metadatamat.urtprctile >= 0.25 & metadatamat.urtprctile < 0.5;
    case 3
        primed_con   = metadatamat.prtprctile >= 0.5  & metadatamat.prtprctile < 0.75;
        unprimed_con = metadatamat.urtprctile >= 0.5  & metadatamat.urtprctile < 0.75;
    case 4
        primed_con   = metadatamat.prtprctile >= 0.75;
        unprimed_con = metadatamat.urtprctile >= 0.75;
end
end
