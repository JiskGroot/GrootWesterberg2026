function [mua_data_corrected, mua_data_saccade_corrected] = baseline_correction(mua_data, mua_data_saccade, baseline_idx, method)
% Compute mean baseline activity per channel and trial
if strcmp(method, 'mean')
    baseline_mean = mean(mua_data(baseline_idx), 'omitnan');
else
    baseline_mean = median(mua_data(baseline_idx), 'omitnan');
end
% Subtract baseline mean from the entire time series
mua_data_corrected = mua_data - baseline_mean;
mua_data_saccade_corrected = mua_data_saccade - baseline_mean;
end