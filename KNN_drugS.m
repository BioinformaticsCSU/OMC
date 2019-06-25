function [y_r] = KNN_drugS(DR_mat, R_mat, K)
%% KNN_drugS:  the function of KNN Preprocessing based on drug similarity matrix
% Usage: [y_r] = KNN_drugS(DR_mat, R_mat, K)
%
% Inputs:
%        DR_mat           - disease-drug association matrix.
%        R_mat            - drug similarity matrix.
%        K                - a positive integer (K nearest neighbor drugs).
%
% Outputs:
%        y_r              - the corresponding column vectors of KNN Preprocessing.

[rows,cols] = size(DR_mat);
y_r = zeros(rows, cols);
col_no = find(sum(DR_mat, 1) == 0);

for i = 1 : length(col_no)
    R_mat(col_no(i), col_no(i)) = 0;
    [sort_d, idx_d] = sort(R_mat(:, col_no(i)), 'descend');
    sum_d = sum(sort_d(1 : K, 1));
    
    if  (sum_d == 0)
        warning('It has no K nearest neighbors!!!')
    end
    
    for j = 1 : K
        y_r(:, col_no(i)) =  y_r(:, col_no(i)) + sort_d(j, 1) * DR_mat(:, idx_d(j, 1));
    end
    y_r(:, col_no(i)) = y_r(:, col_no(i)) / sum_d;
end