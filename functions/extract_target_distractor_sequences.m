function seqMeta = extract_target_distractor_sequences(trialData, session_Num)
% seqMeta = 1 x nOccurences cell array
% each cell: table with 15 * nSeqRows rows (15 channels, consecutive trials),
% and all metadata columns from trialData.

    % Select this session only
    thisSess = trialData(trialData.session == session_Num, :);

    % Sort by trial, then channel to make sure order is consistent
    thisSess = sortrows(thisSess, {'trial','channel'});
    
    trials  = thisSess.trial;
    chans   = thisSess.channel;
    target  = thisSess.target;
    primed  = thisSess.primed;
    
    seqMeta = {};              % will grow: 1 x nOccurences cell
    nRows   = height(thisSess);
    k       = 1;
    
    while k <= nRows
        % look for a target, primed==1
        if target(k) == 1 && primed(k) == 1
            thisTrial = trials(k);
            thisChan  = chans(k);
            
            % require that all 15 channels are present at this trial
            idxTargAllChan = find(trials == thisTrial & chans >= 1 & chans <= 15);
            if numel(idxTargAllChan) ~= 15
                k = k + 1;
                continue
            end
            
            % now build the sequence in trial dimension:
            % we need: trial-1 (all 15 ch), then trial, trial+1, ..., until stop
            seqTrials = thisTrial;         % start with the target trial
            tNext     = thisTrial + 1;
            
            while true
                % rows for next trial and same channel as starting channel,
                % to test the target/primed conditions
                idxNextChan1 = trials == tNext & chans == thisChan;
                if ~any(idxNextChan1)
                    break
                end
                if target(idxNextChan1) ~= 0 || primed(idxNextChan1) ~= 1
                    break
                end
                % valid distractor, primed==1 (for that channel), so extend sequence
                seqTrials(end+1) = tNext; %#ok<AGROW>
                tNext = tNext + 1;
            end
            
            % if we only have the starting target trial and no following distractor,
            % you can decide whether to keep this or skip; here we keep it
            % as a sequence of length 1.
            
            % add the preceding trial (if it exists) to the sequence
            preTrial = thisTrial - 1;
            if any(trials == preTrial)
                seqTrials = [preTrial, seqTrials];
            end
            
            % For each trial in seqTrials, collect all 15 channels
            idxSeq = false(nRows,1);
            for tt = seqTrials
                idxSeq = idxSeq | (trials == tt & chans >= 1 & chans <= 15);
            end
            
            seqTable = thisSess(idxSeq, :);
            
            % sort rows in the sequence by trial then channel
            seqTable = sortrows(seqTable, {'trial','channel'});
            
            % store as one occurrence
            seqMeta{end+1} = seqTable; %#ok<AGROW>
            
            % move k forward: skip over all rows from the last trial in this sequence
            % (for this channel)
            k = find(trials > seqTrials(end) | ...
                     (trials == seqTrials(end) & chans > thisChan), 1, 'first');
            if isempty(k)
                break
            end
        else
            k = k + 1;
        end
    end
end
