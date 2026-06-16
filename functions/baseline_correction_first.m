function [mua_data_array, mua_data_saccade] = baseline_correction_first(mua_data_array, mua_data_saccade, baseline_index, task_data)
% Performs baseline correction on multiunit data for blocks of trials
%
% Inputs:
% - mua_data_array: [channels x time x trials] MUA aligned to array onset
% - mua_data_saccade: [channels x time x trials] MUA aligned to saccade
% - baseline_index: indices for baseline timepoints
% - task_data: table with field .block_trial_count (length = trials)
%
% Output:
% - Baseline-corrected mua_data_array and mua_data_saccade

% Get total number of trials and channels
[channels, ~, num_trials] = size(mua_data_array);
block_trial_count = task_data.block_trial_count;

% Find indices of first trials in each block
first_trial_indices = find(block_trial_count == 1);

% Loop through each block
for i = 1:length(first_trial_indices)
    first_idx = first_trial_indices(i);

    % Determine last trial in the block
    if i < length(first_trial_indices)
        last_idx = first_trial_indices(i+1) - 1;
    else
        last_idx = num_trials;
    end

    % For each channel, calculate baseline from the first trial of the block
    for ch = 1:channels
        baseline = nanmean(mua_data_array(ch, baseline_index, first_idx), 2); % mean over baseline timepoints

        % Subtract baseline from all trials in the block (for both alignments)
        mua_data_array(ch, :, first_idx:last_idx) = ...
            mua_data_array(ch, :, first_idx:last_idx) - baseline;

        mua_data_saccade(ch, :, first_idx:last_idx) = ...
            mua_data_saccade(ch, :, first_idx:last_idx) - baseline;
    end
end

end
