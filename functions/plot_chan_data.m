function fig = plot_chan_data(chan_data,epoch, plotMode, relative, save_bool, save_dir_PNG, save_dir_SVG, avg_compartments)
% chan_data{session, compartment, sequence} = column vector (n x 1), or []
% plotMode: 'all'     -> pool all sessions per compartment & sequence
%           'each'    -> show individual data points per session
%           'session' -> one mean dot per session & sequence (used for Fig 4D)
% relative: 0 -> absolute values
%           1-4 -> subtract that sequence from all others per session
% avg_compartments: 0 -> per-compartment; 1 -> average across compartments first
%
% Saves figures internally when save_bool == 1.

parametric = 0;
anova      = 0;
ttest      = 0;
lines      = 0;

if nargin < 2, plotMode = 'all'; end
if nargin < 3, relative = 0;    end
if nargin < 7, avg_compartments = 0; end

sessions_to_take = 1:29;

[nSessions, nCompartments, nSeq] = size(chan_data);
assert(nCompartments == 3, 'Expected 3 compartments.');
assert(nSeq == 5,          'Expected 5 sequences.');
assert(relative >= 0 && relative <= nSeq, 'relative must be 0..nSeq');

if relative > 0
    chan_data = make_relative(chan_data, nSessions, nCompartments, nSeq, relative);
end

if avg_compartments
    chan_data      = average_compartments(chan_data, nSessions, nCompartments, nSeq);
    nCompartments  = 1;
end

% Colors: first sequence blue, rest dark->light grey
nTrials  = nSeq;
colors   = zeros(nTrials, 3);
colors(1,:) = [0 0 1];
if nTrials > 1
    steps = linspace(0.3, 0.8, nTrials-1)';
    colors(2:end,:) = repmat(steps, 1, 3);
end

if avg_compartments
    compartmentNames = {'AvgCompartments'};
else
    compartmentNames = {'Upper','Middle','Lower'};
end

switch lower(plotMode)

    case 'all'
        for c = 1:nCompartments
            fig = figure; hold on;
            allData = []; allGroup = []; allSess = []; allIdx = [];
            for sess = 1:nSessions
                for s = 1:nSeq
                    x = chan_data{sess, c, s};
                    if isempty(x), continue; end
                    vals = x(:); n = numel(vals);
                    allData  = [allData;  vals];           %#ok<AGROW>
                    allGroup = [allGroup; s*ones(n,1)];    %#ok<AGROW>
                    allSess  = [allSess;  sess*ones(n,1)]; %#ok<AGROW>
                    allIdx   = [allIdx;   (1:n)'];         %#ok<AGROW>
                end
            end
            if ~isempty(allData)
                boxplot(allData, allGroup, 'Positions', 1:nSeq, 'Colors', 'k', 'Symbol', '');
                xJ = nan(size(allData));
                for s = 1:nSeq
                    idx = (allGroup == s); y = allData(idx);
                    x_j = (rand(size(y))-0.5)*0.1;
                    xJ(idx) = s + x_j;
                    scatter(xJ(idx), y, 20, colors(s,:), 'filled', ...
                        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
                end
                for sess = 1:nSessions
                    sessMask = (allSess == sess);
                    if ~any(sessMask), continue; end
                    for k = unique(allIdx(sessMask)).'
                        mask = sessMask & (allIdx == k);
                        if nnz(mask) < 2, continue; end
                        [x_s, ord] = sort(xJ(mask));
                        y_s = allData(mask); y_s = y_s(ord);
                        if lines
                            plot(x_s, y_s, '-', 'Color', [0.7 0.7 0.7 0.2], 'LineWidth', 0.5);
                        end
                    end
                end
                xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
                xlabel('Sequence');
                if relative > 0; ylabel(sprintf('Value (relative to Seq %d)', relative));
                else;            ylabel('Corrected MUA (a.u.)'); end
                title(sprintf('Compartment: %s (all sessions pooled)', compartmentNames{c}));
                box on;
            end
            hold off;
        end

    case 'each'
        comp_loop = 2;
        if avg_compartments, comp_loop = 1; end
        for c = comp_loop
            for sess = 1:nSessions
                if isempty(chan_data{sess, c, 1}), continue; end
                fig = figure; hold on;
                sessionData = []; sessionGroup = [];
                for s = 1:nSeq
                    x = chan_data{sess, c, s};
                    if ~isempty(x)
                        vals = x(:);
                        sessionData  = [sessionData;  vals];              %#ok<AGROW>
                        sessionGroup = [sessionGroup; s*ones(numel(vals),1)]; %#ok<AGROW>
                    end
                end
                if ~isempty(sessionData)
                    boxplot(sessionData, sessionGroup, 'Positions', 1:nSeq, 'Colors', 'k', 'Symbol', '');
                    nT = 0;
                    for s = 1:nSeq
                        x = chan_data{sess, c, s};
                        if ~isempty(x); nT = max(nT, numel(x(:))); break; end
                    end
                    for trial = 1:nT
                        tx = []; ty = [];
                        for s = 1:nSeq
                            x = chan_data{sess, c, s};
                            if ~isempty(x) && trial <= numel(x(:))
                                tx = [tx; s];          %#ok<AGROW>
                                ty = [ty; x(trial)];   %#ok<AGROW>
                            end
                        end
                        if numel(tx) >= 2 && lines
                            plot(tx, ty, '-', 'Color', [0.7 0.7 0.7 0.2], 'LineWidth', 0.7);
                        end
                    end
                    for s = 1:nSeq
                        x = chan_data{sess, c, s};
                        if isempty(x), continue; end
                        y = x(:); xj = (rand(size(y))-0.5)*0.1;
                        scatter(s + xj, y, 20, colors(s,:), 'filled', ...
                            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');
                    end
                    xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
                    xlabel('Sequence');
                    if relative > 0; ylabel(sprintf('Value (relative to Seq %d)', relative));
                    else;            ylabel('Corrected MUA (a.u.)'); end
                    title(sprintf('Compartment: %s - Session %d', compartmentNames{c}, sess));
                    box on;
                end
                hold off;
            end
        end

    case 'session'
        for c = 1:nCompartments
            fig = figure; hold on;
            allData     = [];
            allGroup    = [];
            all_seq_means = nan(nSessions, 1);
            group_means   = nan(nSessions, nSeq);

            for sess = sessions_to_take
                temp = [];
                for s = 1:nSeq
                    temp = [temp; chan_data{sess, c, s}]; %#ok<AGROW>
                end
                all_seq_means(sess) = nanmean(temp);
            end

            for s = 1:nSeq
                for sess = sessions_to_take
                    x = chan_data{sess, c, s};
                    if ~isempty(x)
                        vals = x(:) - all_seq_means(sess);
                        group_means(sess,s) = nanmean(vals);
                        allData  = [allData;  vals];              %#ok<AGROW>
                        allGroup = [allGroup; s*ones(numel(vals),1)]; %#ok<AGROW>
                    end
                end
            end

            if ~isempty(allData)
                if parametric
                    groups_std = nanstd(group_means, 1);
                    errorbar(1:nSeq, nanmean(group_means, 1), groups_std, 'k.', 'LineWidth', 1.5, 'MarkerSize', 8);
                else
                    boxplot(group_means, 'Positions', 1:nSeq, 'Colors', 'k', 'Symbol', '');
                end

                for sess = sessions_to_take
                    x_seq  = 1:nSeq;
                    y_mean = nan(1, nSeq);
                    for s = 1:nSeq
                        x = chan_data{sess, c, s};
                        if ~isempty(x)
                            y_mean(s) = mean(x(:)) - all_seq_means(sess);
                        end
                    end
                    valid = ~isnan(y_mean);
                    if lines && any(valid)
                        plot(x_seq(valid), y_mean(valid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
                    end
                    for s = 1:nSeq
                        if isnan(y_mean(s)), continue; end
                        scatter(s + (rand-0.5)*0.1, y_mean(s), 40, colors(s,:), ...
                            'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'k');
                    end
                end

                if anova
                    Y  = nan(nSessions, nSeq);
                    for sess = sessions_to_take
                        for s = 1:nSeq
                            x = chan_data{sess, c, s};
                            if ~isempty(x); Y(sess,s) = nanmean(x(:)); end
                        end
                    end
                    Yd = Y(:, 2:nSeq);
                    validSess = all(~isnan(Yd), 2);
                    Yd = Yd(validSess, :);
                    [nSubj, nCond] = size(Yd);
                    if nCond > 1 && exist('fitrm', 'file')
                        T  = array2table(Yd, 'VariableNames', compose('D%d', 1:nCond));
                        T.Subj = (1:nSubj)';
                        within = table((2:nSeq)', 'VariableNames', {'Seq'});
                        rm = fitrm(T, sprintf('D1-D%d ~ 1', nCond), 'WithinDesign', within);
                        ranovatbl = ranova(rm, 'WithinModel', 'Seq');
                        p_anova = ranovatbl.pValue(1);
                    else
                        p_anova = NaN;
                    end
                    yl = ylim;
                    xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
                    xlabel('Sequence');
                    if relative > 0; ylabel(sprintf('Value (relative to Seq %d)', relative));
                    else;            ylabel('Corrected MUA (a.u.)'); end
                    if isnan(p_anova)
                        title(sprintf('Compartment: %s (distractor ANOVA n/a)', compartmentNames{c}));
                    else
                        title(sprintf('Compartment: %s (distractor ANOVA p=%.3g)', compartmentNames{c}, p_anova));
                    end
                    ylim([-0.3, 0.3]); yticks([-0.3,-0.15, 0, 0.15, .3]);
                    set(gca, 'TickDir', 'out');

                elseif ttest
                    Y = group_means;
                    pairs = [2 5];
                    pvals = nan(size(pairs,1),1);
                    for k = 1:size(pairs,1)
                        a = pairs(k,1); b = pairs(k,2);
                        valid = ~isnan(Y(:,a)) & ~isnan(Y(:,b));
                        if nnz(valid) >= 2
                            [pvals(k),~,~] = signrank(Y(valid,a), Y(valid,b));
                        end
                    end
                    pvals_corr = pvals * numel(pvals);
                    sig = pvals_corr < 0.05;
                    yl = ylim; y_base = yl(2); y_step = 0.05*range(yl); cur_h = y_base;
                    for k = 1:size(pairs,1)
                        if ~sig(k), continue; end
                        a = pairs(k,1); b = pairs(k,2);
                        xm = (a+b)/2; y_bar = cur_h + y_step;
                        if lines
                            plot([a b], [y_bar y_bar], 'k-', 'LineWidth', 1);
                            plot([a a], [y_bar-0.01*range(yl) y_bar], 'k-', 'LineWidth', 1);
                            plot([b b], [y_bar-0.01*range(yl) y_bar], 'k-', 'LineWidth', 1);
                            text(xm, y_bar+0.01*range(yl), '*', ...
                                'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14);
                        end
                        cur_h = y_bar;
                    end
                    ylim([yl(1), cur_h + 0.1*range(yl)]);
                    xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
                    xlabel('Sequence');
                    if relative > 0; ylabel(sprintf('Value (relative to Seq %d)', relative));
                    else;            ylabel('Corrected MUA (a.u.)'); end
                    title(sprintf('Compartment: %s (session means)', compartmentNames{c}));
                    ylim([-0.4, .4]); yticks([-0.4,-0.2,0,.2,.4]);
                    set(gca, 'TickDir', 'out');

                else
                    % Page's L test on distractor sequences (2..nSeq)
                    valid_rows = ~isnan(group_means(:,1));
                    [L, p_value] = pagesLTest(group_means(valid_rows, 2:end));
                    fprintf('Compartment %s — L = %.2f, p = %.4f\n', compartmentNames{c}, L, p_value);
                end

                xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
                xlabel('Sequence');
                if relative > 0; ylabel(sprintf('Value (relative to Seq %d)', relative));
                else;            ylabel('Corrected MUA (a.u.)'); end
                title(sprintf('Compartment: %s (session means)', compartmentNames{c}));
                ylim([-0.4, .4]); yticks([-0.4,-0.2,0,.2,.4]);
                set(gca, 'TickDir', 'out');
                set(gca, 'Box', 'off');
            end
            hold off;

            if save_bool
                if strcmp(epoch, 'late')
                    if avg_compartments
                        fname_png = 'neural_adap_stats_late_avgComp.png';
                        fname_svg = 'neural_adap_stats_late_avgComp.svg';
                    else
                        fname_png = sprintf('neural_adap_stats_late_%d.png', c);
                        fname_svg = sprintf('neural_adap_stats_late_%d.svg', c);
                    end
                    saveas(fig, fullfile(save_dir_PNG, fname_png));
                    saveas(fig, fullfile(save_dir_SVG, fname_svg));
                else
                    if avg_compartments
                        fname_png = 'neural_adap_stats_early_avgComp.png';
                        fname_svg = 'neural_adap_stats_early_avgComp.svg';
                    else
                        fname_png = sprintf('neural_adap_stats_early_%d.png', c);
                        fname_svg = sprintf('neural_adap_stats_early_%d.svg', c);
                    end
                    saveas(fig, fullfile(save_dir_PNG, fname_png));
                    saveas(fig, fullfile(save_dir_SVG, fname_svg));
                end

            end
        end

    otherwise
        error('Unknown plotMode. Use ''all'', ''each'' or ''session''.');
end
end

% Subtract the specified baseline sequence from all sequences per session/compartment
function chan_data_rel = make_relative(chan_data, nSessions, nCompartments, nSeq, baseline_seq)
chan_data_rel = chan_data;
for sess = 1:nSessions
    for c = 1:nCompartments
        baseline = chan_data{sess, c, baseline_seq};
        if isempty(baseline), continue; end
        baseline = baseline(:); nT = numel(baseline);
        for s = 1:nSeq
            x = chan_data{sess, c, s};
            if isempty(x), continue; end
            x = x(:); nC = min(numel(x), nT);
            chan_data_rel{sess, c, s} = x(1:nC) - baseline(1:nC);
        end
    end
end
end

% Average across compartments, producing a {nSessions, 1, nSeq} cell array
function chan_data_avg = average_compartments(chan_data, nSessions, nCompartments, nSeq)
chan_data_avg = cell(nSessions, 1, nSeq);
for sess = 1:nSessions
    for s = 1:nSeq
        vecs = {};
        for c = 1:nCompartments
            x = chan_data{sess, c, s};
            if ~isempty(x); vecs{end+1} = x(:); end %#ok<AGROW>
        end
        if isempty(vecs); chan_data_avg{sess, 1, s} = []; continue; end
        nT = min(cellfun(@numel, vecs));
        M  = zeros(nT, numel(vecs));
        for k = 1:numel(vecs); M(:,k) = vecs{k}(1:nT); end
        chan_data_avg{sess, 1, s} = nanmean(M, 2);
    end
end
end
