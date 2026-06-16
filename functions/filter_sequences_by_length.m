function allSeqMeta_filt = filter_sequences_by_length(allSeqMeta, nTrials, nOccurences, save_dir_PNG, save_dir_SVG)
% allSeqMeta_filt = filter_sequences_by_length(allSeqMeta, nTrials, nOccurences, save_dir_PNG, save_dir_SVG)
%
% allSeqMeta    : 1xN cell, each cell = 1 x nOccurencesSession cell array of tables
% nTrials       : desired number of task trials in sequence (e.g. 5 for T,D1..D4)
%                 function will require (nTrials + 1) trials including pre-trial
%                 and 15 channels per trial.
% nOccurences   : minimum number of remaining sequences a session must have;
%                 if fewer, that session is set to NaN.
% save_dir_PNG  : directory for PNG output (Figure 4B histogram)
% save_dir_SVG  : directory for SVG output
%
% Output:
% allSeqMeta_filt: same structure as allSeqMeta, but
%                  - sequences with fewer than (nTrials+1)*15 rows removed
%                  - sequences with more rows truncated to exactly (nTrials+1)*15 rows
%                  - sessions with < nOccurences sequences replaced by NaN.
%                  Also saves Figure 4B (session-count bar chart) internally.

    minRows   = (nTrials + 1) * 15;
    nSessions = numel(allSeqMeta);
    allSeqMeta_filt     = cell(size(allSeqMeta));
    total_sequences     = 0;
    total_sessions      = 0;
    seqCountsPerSession = zeros(1, nSessions);

    for s = 1:nSessions
        seqs = allSeqMeta{s};
        if isempty(seqs)
            allSeqMeta_filt{s}      = NaN;
            seqCountsPerSession(s)  = 0;
            continue
        end
        seqHeights = cellfun(@height, seqs);
        keepMask   = seqHeights >= minRows;
        seqs_trunc = cellfun(@(t) t(1:minRows, :), seqs(keepMask), 'UniformOutput', false);
        nRemain    = numel(seqs_trunc);
        seqCountsPerSession(s) = nRemain;
        if nRemain < nOccurences
            allSeqMeta_filt{s} = NaN;
        else
            allSeqMeta_filt{s}  = seqs_trunc;
            total_sequences     = total_sequences + nRemain;
            total_sessions      = total_sessions  + 1;
        end
    end
    sprintf("Total sequences: %d, over %d sessions", total_sequences, total_sessions);

    % ---- Figure 4B: session-count bar chart with monkey color-coding ----
    nCa = min(21, nSessions);                     % sessions 1-21 = Monkey Ca
    nHe = max(0,  nSessions - 21);                % sessions 22-end = Monkey He

    caColor      = [0.2 0.4 0.8];
    heColor      = [0.8 0.4 0.2];
    defaultColor = [0.6 0.6 0.6];

    fig = figure; hold on;
    b = bar(1:nSessions, seqCountsPerSession, 'FaceColor', 'flat');
    b.CData = repmat(defaultColor, nSessions, 1);
    if nCa > 0
        b.CData(1:nCa, :) = repmat(caColor, nCa, 1);
    end
    if nHe > 0
        b.CData(22:min(21+nHe, nSessions), :) = repmat(heColor, nHe, 1);
    end

    yline(nOccurences, 'r--', 'LineWidth', 1.5);
    xlabel('Session');
    ylabel('Number of occurrences of full sequence');
    title('Full-sequence occurrences per session');
    ylim([0 75]);   yticks([0, 25, 50, 75]);
    xlim([0, nSessions + 1]);
    xticks(1:nSessions);
    set(gca, 'TickDir', 'out');

    % Legend with totals
    totalCa = sum(seqCountsPerSession(1:nCa));
    totalHe = sum(seqCountsPerSession(22:min(21+nHe, nSessions)));

    legendEntriesValid = {};
    legendColorsValid  = [];
    if nCa > 0
        legendEntriesValid{end+1} = sprintf('Monkey Ca (n = %d)', totalCa);
        legendColorsValid = [legendColorsValid; caColor];
    end
    if nHe > 0
        legendEntriesValid{end+1} = sprintf('Monkey He (n = %d)', totalHe);
        legendColorsValid = [legendColorsValid; heColor];
    end

    if ~isempty(legendEntriesValid)
        legHandles = gobjects(size(legendColorsValid,1), 1);
        for i = 1:numel(legHandles)
            legHandles(i) = plot(NaN, NaN, 's', ...
                'MarkerFaceColor', legendColorsValid(i,:), ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1.5);
        end
        legend(legHandles, legendEntriesValid, 'Location', 'best', 'FontSize', 9);
    end

    hold off;
    saveas(fig, fullfile(save_dir_PNG, 'neural_adaptation_session_count.png'));
    saveas(fig, fullfile(save_dir_SVG, 'neural_adaptation_session_count.svg'));
end
