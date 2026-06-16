function [p] = running_rs_null(data_mat1, data_bsl,alpha)
p = nan(1, size(data_mat1,2));
for i = 1:size(data_mat1, 2)
    data_1 = data_mat1(:,i);
    if any(~isnan(data_1))
        p(i) = ranksum(data_1, data_bsl, 'alpha', alpha);
    end
end

end