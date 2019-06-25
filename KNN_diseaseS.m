function [y_d] = KNN_diseaseS(DR_mat, D_mat, K)
%% KNN_diseaseS:  the function of KNN Preprocessing based on disease similarity matrix
% Usage: [y_d] = KNN_diseaseS(DR_mat, D_mat, K)
%
% Inputs:
%        DR_mat           - disease-drug association matrix.
%        D_mat            - disease similarity matrix.
%        K                - a positive integer (K nearest neighbor diseases).
%
% Outputs:
%        y_d              - the corresponding row vectors of KNN Preprocessing.

[rows,cols] = size(DR_mat);
y_d = zeros(rows, cols);
row_no = find(sum(DR_mat, 2) == 0);

for i = 1 : length(row_no)
    D_mat(row_no(i), row_no(i)) = 0;
    [sort_d, idx_d] = sort(D_mat(row_no(i), :), 'descend');
    sum_d = sum(sort_d(1, 1 : K));
    
    if  (sum_d == 0)
        warning('It has no K nearest neighbors!!!')
    end
    
    for j = 1 : K
        y_d(row_no(i), :) =  y_d(row_no(i), :)+ sort_d(1, j) * DR_mat(idx_d(1, j), :);
    end
    y_d(row_no(i), :) = y_d(row_no(i), :) / sum_d;
end