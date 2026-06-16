function [L, p_value] = pagesLTest(group_means)
% Page's L test for a monotonically decreasing trend across conditions.
%
% group_means : n x k matrix (n observations, k ordered conditions)
% L           : Page's L statistic
% p_value     : one-sided p-value (normal approximation)

[n, k] = size(group_means);

% Rank each row (handles ties)
ranks = zeros(size(group_means));
for i = 1:n
    ranks(i,:) = tiedrank(group_means(i,:));
end

% Sum ranks per condition column
R = sum(ranks, 1);

% Decreasing weights (k, k-1, ..., 1)
w = k:-1:1;

% Page's L statistic
L = sum(w .* R);

% Normal approximation
mu_L    = n * k * (k+1)^2 / 4;
sigma_L = sqrt(n * k^2 * (k+1) * (k^2-1) / 144);
z       = (L - mu_L) / sigma_L;

% One-sided (decreasing trend)
p_value = 1 - normcdf(z);
end
