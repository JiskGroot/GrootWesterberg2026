function [MUA_datamat, MUA_metadatamat, CSD_datamat, CSD_metadatamat, params] = preprocess_neural_data(params, varargin)
% PREPROCESS_NEURAL_DATA: Optimized single-step preprocessing pipeline
% Loads raw data and outputs final datamat/metadatamat structure for both MUA and CSD
%
% Inputs:
%   params - struct with fields:
%       .disk_dir, .monkey_name, .time_data, .cut_off, .trial_threshold
%       .zscore_session_wise, .rel_time_point ('ao' or 'so')
%   varargin - optional name-value pairs:
%       'baseline_correction' - 'norm' (default) or 'first'
%       'z_score' - 1 (default) or 0
%       'smooth_data' - 1 (default) or 0
%       'smooth_strength' - 15 (default)
%
% Outputs:
%   MUA_datamat - [observations x timepoints] matrix
%   MUA_metadatamat - table with trial metadata
%   CSD_datamat - [observations x timepoints] matrix
%   CSD_metadatamat - table with trial metadata (same as MUA)
%   params - updated params struct

%% Parse input arguments
baseline_correction = 'norm';
z_score = 1;
smooth_data = 1;
smooth_strength_MUA = 15;
smooth_strength_CSD = 10;

varStrInd = find(cellfun(@ischar, varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case 'baseline_correction'
            baseline_correction = varargin{varStrInd(iv) + 1};
        case 'z_score'
            z_score = varargin{varStrInd(iv) + 1};
        case 'smooth_data'
            smooth_data = varargin{varStrInd(iv) + 1};
        case 'smooth_strength_MUA'
            smooth_strength_MUA = varargin{varStrInd(iv) + 1};
        case 'smooth_strength_CSD'
            smooth_strength_CSD = varargin{varStrInd(iv) + 1};
    end
end

%% Initialize parameters
disk_dir = params.disk_dir;
monkey_name = params.monkey_name;
baseline_period = [-250, 0];

% Get list of datasets
if strcmp(monkey_name, "both")
    file_list = dir(fullfile(disk_dir, '*-npy'));
elseif strcmp(monkey_name, "C")
    file_list = dir(fullfile(disk_dir, 'C*-npy'));
elseif strcmp(monkey_name, "H")
    file_list = dir(fullfile(disk_dir, 'H*-npy'));
end
dataset_names = {file_list.name};

%% Initialize final output matrices
MUA_datamat = [];
CSD_datamat = [];
metadatamat = table();
obs_ctr = 0;

%% Main loop: Process each session
for session = 1:length(dataset_names)
    dataset_name = dataset_names{session};
    dataset_path = fullfile(disk_dir, dataset_name);

    fprintf('Processing session %d/%d: %s\n', session, length(dataset_names), dataset_name);

    %% Load raw data
    try
        recording_data = matread_csv(fullfile(dataset_path, "recordinginfo.csv"));
        task_data = matread_csv(fullfile(dataset_path, "task.csv"));
        time_data = matread_npy(fullfile(dataset_path, "time_array_ms.npy"));
        behavior_data = matread_npy(fullfile(dataset_path, "behavior.npy"));

        % Load MUA and CSD data based on alignment
        if strcmp(params.rel_time_point, "ao")
            mua_data = matread_npy(fullfile(dataset_path, "mua_array_uv.npy"));
            csd_data = matread_npy(fullfile(dataset_path, "csd_array_naPmm3.npy"));
        else % saccade onset
            mua_data = matread_npy(fullfile(dataset_path, "mua_saccade_uv.npy"));
            csd_data = matread_npy(fullfile(dataset_path, "csd_saccade_naPmm3.npy"));
        end
    catch
        warning('Failed to load data for %s, skipping...', dataset_name);
        continue;
    end

    %% Validate data dimensions
    if size(mua_data, 1) ~= 15 || ndims(mua_data) ~= 3
        fprintf('Skipping %s: incorrect dimensions\n', dataset_name);
        continue;
    end

    % Fix CSD dimensions (remove NaN row)
    csd_data = csd_data(2:16, :, :);

    %% Define baseline and time indices
    baseline_index = (time_data >= baseline_period(1)) & (time_data < baseline_period(2));

    %% Apply baseline correction to first trial of each block (if 'first' method)
    if strcmp(baseline_correction, 'first')
        [mua_data, ~] = baseline_correction_first(mua_data, mua_data, baseline_index, task_data);
        [csd_data, ~] = baseline_correction_first(csd_data, csd_data, time_data >= -50 & time_data < 0, task_data);
    end

    %% Calculate session-wise baseline SD (if using session-wise z-score)
    if params.zscore_session_wise == 1
        bsl_array_mua = mua_data(:, baseline_index, :);
        bsl_array_sd_mua = std(bsl_array_mua, 0, [2,3], 'omitnan');

        bsl_array_csd = csd_data(:, time_data >= -50 & time_data < 0, :);
        bsl_array_sd_csd = std(bsl_array_csd, 0, [2,3], 'omitnan');
    end

    %% Identify target and distractor trials
    trials_target = (task_data.array_target_position_deg == recording_data.receptive_field_position_deg) & ...
        (behavior_data(:, 2) == 1) & ~isnan(behavior_data(:, 1));
    trials_distractor = (task_data.array_target_position_deg ~= recording_data.receptive_field_position_deg) & ...
        (task_data.catch_trial_logical == 0) & (behavior_data(:, 2) == 1) & ~isnan(behavior_data(:, 1));

    %% Process trials by block and target type
    block_counts = unique(task_data.block_trial_count);

    temp_mua_datamat = [];
    temp_csd_datamat = [];
    temp_metadatamat = table();

    for block_num = 1:min(params.cut_off, max(block_counts))
        block_mask = (task_data.block_trial_count == block_num);

        % Process both target (kk=3) and distractor (kk=2) conditions
        for kk = 2:3
            if kk == 2
                trial_mask = trials_distractor & block_mask;
                is_target = 0;
            else
                trial_mask = trials_target & block_mask;
                is_target = 1;
            end

            % Skip if insufficient trials
            ntrial = sum(trial_mask);
            if ntrial < 2
                continue;
            end

            %% Process MUA and CSD data for these trials
            mua_block = mua_data(:, :, trial_mask);
            csd_block = csd_data(:, :, trial_mask);
            rt_block = behavior_data(trial_mask, 1);

            if kk == 2
                color_block = task_data.distractor_color_string(trial_mask);
            else
                color_block = task_data.target_color_string(trial_mask);
            end

            %% Apply preprocessing to MUA
            if params.zscore_session_wise == 1
                mua_processed = process_trials_session_wise(mua_block, rt_block, baseline_index, ...
                    time_data, smooth_data, z_score, smooth_strength_MUA, baseline_correction, ...
                    bsl_array_sd_mua, params.rel_time_point);
            else
                mua_processed = process_trials_trial_wise(mua_block, rt_block, baseline_index, ...
                    time_data, smooth_data, z_score, smooth_strength_MUA, baseline_correction, params.rel_time_point);
            end

            %% Apply preprocessing to CSD
            csd_baseline_index = time_data >= -50 & time_data < 0;
            if params.zscore_session_wise == 1
                csd_processed = process_trials_session_wise(csd_block, rt_block, csd_baseline_index, ...
                    time_data, 1, 0, smooth_strength_CSD, baseline_correction, bsl_array_sd_csd, params.rel_time_point);
            else
                csd_processed = process_trials_trial_wise(csd_block, rt_block, csd_baseline_index, ...
                    time_data, 1, 0, smooth_strength_CSD, baseline_correction, params.rel_time_point);
            end

            %% Concatenate processed data across trials (3rd dimension)
            temp_mua_datamat = cat(3, temp_mua_datamat, mua_processed);
            temp_csd_datamat = cat(3, temp_csd_datamat, csd_processed);

            %% Create metadata for these trials
            temp_meta = table();
            temp_meta.target = zeros(ntrial, 1) + is_target;
            temp_meta.rt = rt_block;
            temp_meta.color = color_block;
            temp_meta.primed = zeros(ntrial, 1) + double(block_num > 2);
            temp_meta.sincechange = zeros(ntrial, 1) + block_num;
            temp_meta.channel = ones(ntrial, 1); % Will be replicated later
            temp_meta.session = zeros(ntrial, 1) + session;
            temp_meta.trial = find(trial_mask);
            temp_meta.obs = (obs_ctr + 1 : obs_ctr + ntrial)';

            obs_ctr = obs_ctr + ntrial;
            temp_metadatamat = [temp_metadatamat; temp_meta]; %#ok<AGROW>
        end
    end

    %% Skip session if insufficient trials
    if isempty(temp_metadatamat) || height(temp_metadatamat) < params.trial_threshold
        fprintf('Skipping session %s: insufficient trials\n', dataset_name);
        continue;
    end

    %% Calculate RT percentiles within session
    [~, p] = sort(temp_metadatamat.rt, 'descend');
    r = 1:length(temp_metadatamat.rt);
    temp_metadatamat.rtprctile = zeros(length(r), 1);
    temp_metadatamat.rtprctile(p) = rescale(r)';

    % Unprimed RT percentiles
    temp_metadatamat.urtprctile = nan(size(temp_metadatamat.rtprctile));
    temp_urt = temp_metadatamat.rt(temp_metadatamat.primed == 0);
    [~, p] = sort(temp_urt, 'descend');
    r = 1:numel(temp_urt);
    temp_urt_rescaled = zeros(size(temp_urt));
    temp_urt_rescaled(p) = rescale(r)';
    temp_metadatamat.urtprctile(temp_metadatamat.primed == 0) = temp_urt_rescaled;

    % Primed RT percentiles
    temp_metadatamat.prtprctile = nan(size(temp_metadatamat.rtprctile));
    temp_prt = temp_metadatamat.rt(temp_metadatamat.primed == 1);
    if ~isempty(temp_prt)
        [~, p] = sort(temp_prt, 'descend');
        r = 1:numel(temp_prt);
        temp_prt_rescaled = zeros(size(temp_prt));
        temp_prt_rescaled(p) = rescale(r)';
        temp_metadatamat.prtprctile(temp_metadatamat.primed == 1) = temp_prt_rescaled;
    end

    %% Replicate metadata across all 15 channels
    n_meta    = height(temp_metadatamat);
    full_meta = repmat(temp_metadatamat, 15, 1);
    full_meta.channel = repelem((1:15)', n_meta);

    %% Reshape data: from [channels x time x trials] to [observations x time]
    % observations = channels * trials, channel-major order
    n_channels     = size(temp_mua_datamat, 1);
    n_time_local   = size(temp_mua_datamat, 2);
    n_trials_local = size(temp_mua_datamat, 3);
    session_mua_data = reshape(permute(temp_mua_datamat, [3 1 2]), n_channels * n_trials_local, n_time_local);
    session_csd_data = reshape(permute(temp_csd_datamat, [3 1 2]), n_channels * n_trials_local, n_time_local);

    %% Add color preference
    if params.col_analysis
        full_meta = feature_sel(session_mua_data, full_meta, params);
    end
    
    %% Append to final matrices
    MUA_datamat = [MUA_datamat; session_mua_data]; %#ok<AGROW>
    CSD_datamat = [CSD_datamat; session_csd_data]; %#ok<AGROW>
    metadatamat = [metadatamat; full_meta]; %#ok<AGROW>
end

%% Return outputs
MUA_metadatamat = metadatamat;
CSD_metadatamat = metadatamat; % Same metadata for both

fprintf('\nPreprocessing complete: %d observations across %d sessions\n', height(metadatamat), session);

end

%% ========== HELPER FUNCTIONS ==========

function data_processed = process_trials_trial_wise(data, rt, baseline_index, time_data, ...
    do_smooth, do_zscore, smooth_strength, baseline_correction, rel_time_point)
% Process trials with trial-wise baseline correction and z-scoring

[n_channels, ~, n_trials] = size(data);
data_processed = nan(size(data));

if do_smooth
    b = ones(1, smooth_strength) / smooth_strength;
    a = 1;
end

for trial_i = 1:n_trials
    response_time = rt(trial_i);
    clip_index = (time_data >= (response_time - 10));

    for channel_i = 1:n_channels
        data_temp = squeeze(data(channel_i, :, trial_i));

        % Skip broken channels
        if all(data_temp == data_temp(1))
            continue;
        end

        % Smoothing
        if do_smooth
            data_temp = filtfilt(b, a, data_temp);
        end

        % Z-score (trial-wise)
        if do_zscore
            baseline_std = std(data_temp(baseline_index), 'omitnan');
            if ~isnan(baseline_std) && baseline_std > 0
                data_temp = data_temp / baseline_std;
            end
        end

        % Baseline correction (subtract mean)
        if strcmp(baseline_correction, 'norm')
            baseline_mean = mean(data_temp(baseline_index), 'omitnan');
            data_temp = data_temp - baseline_mean;
        end

        % Clip response artifact for array onset
        if strcmp(rel_time_point, "ao")
            data_temp(clip_index) = NaN;
        end

        data_processed(channel_i, :, trial_i) = data_temp;
    end
end
end

function data_processed = process_trials_session_wise(data, rt, baseline_index, time_data, ...
    do_smooth, do_zscore, smooth_strength, baseline_correction, bsl_array_sd, rel_time_point)
% Process trials with session-wise z-scoring

[n_channels, ~, n_trials] = size(data);
data_processed = nan(size(data));

if do_smooth
    b = ones(1, smooth_strength) / smooth_strength;
    a = 1;
end

for trial_i = 1:n_trials
    response_time = rt(trial_i);
    clip_index = (time_data >= (response_time - 10));

    for channel_i = 1:n_channels
        data_temp = squeeze(data(channel_i, :, trial_i));

        % Skip broken channels
        if all(data_temp == data_temp(1))
            continue;
        end

        % Z-score (session-wise)
        if do_zscore
            if ~isnan(bsl_array_sd(channel_i)) && bsl_array_sd(channel_i) > 0
                data_temp = data_temp / bsl_array_sd(channel_i);
            end
        end

        % Smoothing
        if do_smooth
            data_temp = filtfilt(b, a, data_temp);
        end

        % Baseline correction (subtract mean)
        if strcmp(baseline_correction, 'norm')
            baseline_mean = mean(data_temp(baseline_index), 'omitnan');
            data_temp = data_temp - baseline_mean;
        end

        % Clip response artifact for array onset
        if strcmp(rel_time_point, "ao")
            data_temp(clip_index) = NaN;
        end

        data_processed(channel_i, :, trial_i) = data_temp;
    end
end
end

function [data_array, data_saccade] = baseline_correction_first(data_array, data_saccade, baseline_index, task_data)
% Baseline correction using first trial of each block

[channels, ~, num_trials] = size(data_array);
block_trial_count = task_data.block_trial_count;
first_trial_indices = find(block_trial_count == 1);

for i = 1:length(first_trial_indices)
    first_idx = first_trial_indices(i);

    if i < length(first_trial_indices)
        last_idx = first_trial_indices(i+1) - 1;
    else
        last_idx = num_trials;
    end

    for ch = 1:channels
        baseline = mean(data_array(ch, baseline_index, first_idx), 2, 'omitnan');
        data_array(ch, :, first_idx:last_idx) = data_array(ch, :, first_idx:last_idx) - baseline;
        data_saccade(ch, :, first_idx:last_idx) = data_saccade(ch, :, first_idx:last_idx) - baseline;
    end
end
end
function full_meta = feature_sel(session_mua_data, full_meta, params)
% Allocate time of interest
time_data = params.time_data;
time_data_window = time_data >= 60 & time_data <= 160; % from nat comms
n_permutations = 1000; % 1000 permutions
n_units = 100; % of 100 units each

% Create new columns
full_meta.pref = nan(size(full_meta,1),1);
full_meta.pref_pr = nan(size(full_meta,1),1);

% Find red and green trials
red_trials = unique(full_meta.trial(strcmp(full_meta.color, 'red'))); % find red trials (unique because 15 instances - 1 per channel)
green_trials = unique(full_meta.trial(strcmp(full_meta.color, 'green')));


if length(red_trials) >= n_units && length(green_trials) >= n_units
    % Run population reliability test
    red_count = 0;
    for i_perm = 1:n_permutations
        % draw random trials
        rs_red = red_trials(randperm(length(red_trials), n_units));
        rs_green = green_trials(randperm(length(green_trials), n_units));
        % Sum draws of current permutation
        red_draw = sum(nanmean(session_mua_data(ismember(full_meta.trial, rs_red), time_data_window), 2));
        green_draw = sum(nanmean(session_mua_data(ismember(full_meta.trial, rs_green), time_data_window), 2));
        % Count times red > green
        if red_draw > green_draw
            red_count = red_count + 1;
        end
    end
    % Estimate prob
    prob_red = ((red_count / n_permutations) * 100) - 50; % tuned to red (50), green (-50)
    if prob_red < -2
        full_meta.pref(ismember(full_meta.trial, green_trials)) = 1;
        full_meta.pref(ismember(full_meta.trial, red_trials)) = -1;
    elseif prob_red > 2
        full_meta.pref(ismember(full_meta.trial, green_trials)) = -1;
        full_meta.pref(ismember(full_meta.trial, red_trials)) = 1;
    else
        full_meta.pref = zeros(size(full_meta,1),1);
    end

    % Save selectivity value
    full_meta.pref_pr = repmat(prob_red, size(full_meta,1),1);
end
end