function [mean_out, lower_out, upper_out] = confidence_interval(data_in, methd, lvl)

if nargin < 2 || isempty(methd); methd = 'mean'; end
if nargin < 3 || isempty(lvl);   lvl   = 1.96;   end

if strcmp(methd, 'mean')
    avg_fn = @(x, dim) nanmean(x, dim);
else
    avg_fn = @(x, dim) nanmedian(x, dim);
end

if ndims(data_in) == 3
    % data_in: channel x time x trial
    mean_out  = avg_fn(data_in, 3);
    se        = nanstd(data_in, [], 3) ./ sqrt(size(data_in, 3));
else
    % data_in: observations x time
    mean_out  = avg_fn(data_in, 1);
    se        = nanstd(data_in) ./ sqrt(size(data_in, 1));
end

upper_out = mean_out + lvl * se;
lower_out = mean_out - lvl * se;
end
