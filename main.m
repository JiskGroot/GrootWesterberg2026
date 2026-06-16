% MAIN.M — Data loading and preprocessing
% Run this script first to load and preprocess MUA and CSD data.
% After this, run analysis_MUA.m to reproduce all figures in the paper.
% Loading data (main.m), preprocessing (main.m), and analyses (analysis.m)
% takes approximately 30 min. in total
% Requires: data in NPY/CSV format (PrimateV4cxcol_PrimingOfPopout-npy)

close all

%% Configure paths
disk_dir        = "C:/path/to/datafolder/PrimateV4cxcol_PrimingOfPopout-npy/";
figure_save_dir = "C:/path/to/figure/save/folder";

%% Create output directories
dirs_to_create = {
    fullfile(figure_save_dir, 'MUA', 'PNG')
    fullfile(figure_save_dir, 'MUA', 'SVG')
    fullfile(figure_save_dir, 'MUA', 'PNG', 'lam_trace')
    fullfile(figure_save_dir, 'MUA', 'SVG', 'lam_trace')
    fullfile(figure_save_dir, 'CSD', 'PNG')
    fullfile(figure_save_dir, 'CSD', 'SVG')
};
for i = 1:numel(dirs_to_create)
    if ~exist(dirs_to_create{i}, 'dir')
        mkdir(dirs_to_create{i});
    end
end

%% Add functions folder to path
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Set preprocessing parameters
params.disk_dir             = disk_dir;
params.monkey_name          = "both"; % monkey to use
params.cut_off              = 15; %(15)trial after feature change cut-off 
params.trial_threshold      = 0; % (0)
params.zscore_session_wise  = 0; % (0) 0 = z-score trial wise | 1 = session wise 
params.rel_time_point       = "ao";   % "ao" = array onset | "so" = saccade onset
params.time_data            = matread_npy(fullfile(disk_dir, "C190127-npy/time_array_ms.npy")); % get time from first session
params.figure_save_dir      = figure_save_dir;
params.col_analysis         = 1; % (1) 1 = estimate session-wise color pref. | 0 = not

%% Preprocess data
[MUA_datamat, metadatamat, CSD_datamat, ~, params] = ...
    preprocess_neural_data(params, ...
        'baseline_correction', 'norm', ... % ('norm') how to perform bsl corr (norm | first)
        'z_score',             1,  ... % (1) to z-score data
        'smooth_data',       1,  ... % (1) to smooth data
        'smooth_strength_MUA',     15, ... % (15) movmean win size MUA
        'smooth_strength_CSD',     1); % (1) movmean win size CSD

fprintf('\nPreprocessing complete. Ready to run analyses.\n');

