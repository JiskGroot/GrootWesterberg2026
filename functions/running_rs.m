function [p] = running_rs(data_mat1, data_mat2,alpha)
p = nan(1, size(data_mat1,2));
for i = 1:size(data_mat1, 2)
    data_1 = data_mat1(:,i);
    data_2 = data_mat2(:,i);
    if any(~isnan(data_1)) && any(~isnan(data_2))
        p(i) = ranksum(data_1, data_2, 'alpha', alpha);
    end
end

end