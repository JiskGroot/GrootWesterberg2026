function fig = plot_chan_diff_session(chan_data1, chan_data2, relative, save_bool, save_dir_PNG, save_dir_SVG)
% Session-demean each dataset, compute their difference, and plot per compartment.
%
% Steps:
%   1) For each session+compartment, subtract the mean across all sequences.
%   2) Subtract demeaned dataset2 from demeaned dataset1.
%   3) Plot in 'session' mode: one mean dot per session & sequence.
%
% save_bool = true saves one figure per compartment to save_dir_PNG/SVG.

if nargin < 3, relative  = 0;     end
if nargin < 4, save_bool = false;  end

[nSessions1, nCompartments1, nSeq1] = size(chan_data1);
[nSessions2, nCompartments2, nSeq2] = size(chan_data2);
assert(isequal([nSessions1 nCompartments1 nSeq1], [nSessions2 nCompartments2 nSeq2]), ...
    'chan_data1 and chan_data2 must have the same size.');

nSessions     = nSessions1;
nCompartments = nCompartments1;
nSeq          = nSeq1;

assert(nCompartments == 3, 'Expected 3 compartments.');
assert(relative >= 0 && relative <= nSeq, 'relative must be between 0 and nSeq.');

sessions_to_take = 1:nSessions;
lines      = 0;
parametric = 0;

% ---- 1) Session-demean each dataset ----
demeaned1 = chan_data1;
demeaned2 = chan_data2;

for sess = 1:nSessions
    for c = 1:nCompartments
        all1 = [];
        for s = 1:nSeq
            x1 = chan_data1{sess, c, s};
            if ~isempty(x1); all1 = [all1; x1(:)]; end %#ok<AGROW>
        end
        mu1 = nanmean(all1);

        all2 = [];
        for s = 1:nSeq
            x2 = chan_data2{sess, c, s};
            if ~isempty(x2); all2 = [all2; x2(:)]; end %#ok<AGROW>
        end
        mu2 = nanmean(all2);

        for s = 1:nSeq
            x1 = chan_data1{sess, c, s};
            demeaned1{sess, c, s} = deal_demean(x1, mu1);
            x2 = chan_data2{sess, c, s};
            demeaned2{sess, c, s} = deal_demean(x2, mu2);
        end
    end
end

% ---- 2) Difference of demeaned datasets ----
diff_data = demeaned1;
for sess = 1:nSessions
    for c = 1:nCompartments
        for s = 1:nSeq
            x1 = demeaned1{sess, c, s};
            x2 = demeaned2{sess, c, s};
            if isempty(x1) || isempty(x2)
                diff_data{sess, c, s} = [];
            else
                x1 = x1(:); x2 = x2(:);
                n  = min(numel(x1), numel(x2));
                diff_data{sess, c, s} = x1(1:n) - x2(1:n);
            end
        end
    end
end

% ---- 3) Optional relative transform ----
if relative > 0
    diff_data = make_relative(diff_data, nSessions, nCompartments, nSeq, relative);
end

% ---- Plot ----
colors = zeros(nSeq, 3);
colors(1,:) = [0 0 1];
if nSeq > 1
    steps = linspace(0.3, 0.8, nSeq-1)';
    colors(2:end,:) = repmat(steps, 1, 3);
end
compartmentNames = {'Upper','Middle','Lower'};

for c = 1:nCompartments
    fig = figure; hold on;

    allData     = [];
    allGroup    = [];
    group_means = nan(nSessions, nSeq);

    for s = 1:nSeq
        for sess = sessions_to_take
            x = diff_data{sess, c, s};
            if ~isempty(x)
                vals = x(:);
                group_means(sess,s) = nanmean(vals);
                allData  = [allData;  vals];              %#ok<AGROW>
                allGroup = [allGroup; s*ones(numel(vals),1)]; %#ok<AGROW>
            end
        end
    end

    if ~isempty(allData)
        if parametric
            groups_std = nanstd(group_means, 1);
            errorbar(1:nSeq, nanmean(group_means,1), groups_std, 'k.', 'LineWidth', 1.5, 'MarkerSize', 8);
        else
            boxplot(group_means, 'Positions', 1:nSeq, 'Colors', 'k', 'Symbol', '');
        end

        for sess = sessions_to_take
            x_seq  = 1:nSeq;
            y_mean = group_means(sess,:);
            valid  = ~isnan(y_mean);
            if lines && any(valid)
                plot(x_seq(valid), y_mean(valid), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
            end
            for s = 1:nSeq
                if isnan(y_mean(s)), continue; end
                scatter(s + (rand-0.5)*0.1, y_mean(s), 40, colors(s,:), ...
                    'filled', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'k');
            end
        end

        % Signed-rank test on pair [2 5]
        Y     = group_means;
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
        yl = ylim; cur_h = yl(2); y_step = 0.05*range(yl);
        for k = 1:size(pairs,1)
            if ~sig(k), continue; end
            a = pairs(k,1); b = pairs(k,2);
            xm = (a+b)/2; y_bar = cur_h + y_step;
            if lines
                plot([a b], [y_bar y_bar], 'k-', 'LineWidth', 1);
                plot([a a], [y_bar-0.01*range(yl) y_bar], 'k-', 'LineWidth', 1);
                plot([b b], [y_bar-0.01*range(yl) y_bar], 'k-', 'LineWidth', 1);
                text(xm, y_bar+0.01*range(yl), '*', 'HorizontalAlignment','center', ...
                    'VerticalAlignment','bottom', 'FontSize', 14);
            end
            cur_h = y_bar;
        end
        ylim([yl(1), cur_h + 0.1*range(yl)]);
    end

    xlim([0.5 nSeq+0.5]); xticks(1:nSeq);
    xlabel('Sequence');
    if relative > 0
        ylabel(sprintf('Demeaned diff (relative to Seq %d)', relative));
    else
        ylabel('Session-demeaned corrected MUA diff (a.u.)');
    end
    title(sprintf('Compartment: %s  (data1−sessmean1) − (data2−sessmean2)', compartmentNames{c}));
    ylim([-0.4, .4]); yticks([-0.4,-0.2,0,.2,.4]);
    set(gca, 'TickDir', 'out');
    hold off;

    if save_bool
        saveas(fig, fullfile(save_dir_PNG, sprintf('neural_adap_diff_session_%d.png', c)));
        saveas(fig, fullfile(save_dir_SVG, sprintf('neural_adap_diff_session_%d.svg', c)));
    end
end
end

function out = deal_demean(x, mu)
if isempty(x) || isnan(mu); out = []; else; out = x(:) - mu; end
end

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
