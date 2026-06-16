function p = running_kwt(data_cell)

n_time = size(data_cell{1}, 2);
n_rows = sum(cellfun(@(x) size(x, 1), data_cell));

data  = nan(n_rows, n_time);
group = zeros(n_rows, 1);
row   = 1;
for ii = 1:numel(data_cell)
    n = size(data_cell{ii}, 1);
    data(row:row+n-1, :) = data_cell{ii};
    group(row:row+n-1)   = ii;
    row = row + n;
end

p = nan(1, n_time);
for ii = 1:n_time
    p(ii) = kruskalwallis(data(:,ii), group, "off");
end

end
