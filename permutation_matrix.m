function [ Pmat, auxi] = permutation_matrix( zmatrix, matrix1red, matrix2red )
%

auxi = zeros(1,max(size(zmatrix)));
Pmat = zeros(size(matrix1red,1),size(matrix2red,1));
% construct permutation matrix, 'Pmat'.
for i=1:max(size(zmatrix))
    [~,I] = min(zmatrix(:));
    [I_row, I_col] = ind2sub(size(zmatrix),I);
    Pmat(I_row,I_col) = 1;
    % ¡®auxi¡¯ store the corresponding matching indices
    auxi(I_col) = I_row;
    zmatrix(:,I_col) = 255*ones(size(zmatrix,1),1);
    zmatrix(I_row,:) = 255*ones(1,size(zmatrix,2));
end

end

