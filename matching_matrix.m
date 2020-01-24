function [ zmatrix,matrix1histmod,matrix2histmod ] = matching_matrix(K, matrix1red, matrix2red, helpyn, valormin )
%% Find Zmatrix, using the theory of abs(U)abs(U2) = I or not


[ matrix1histmod,matrix2histmod ] = pair_eigenvector(K, matrix1red, matrix2red, helpyn, valormin );

% ¡®zmatrix¡¯ is the matching matrix in paper in step 4, which compute the difference between eigenvector
zmatrix = 255*ones(size(matrix1histmod,1),size(matrix2histmod,1));
for i=1:size(matrix1histmod,1)
    V = matrix1histmod(i,:);
    for j=1:size(matrix2histmod,1)
        V2 = matrix2histmod(j,:);
        % ¡®zmatrix¡¯
        zmatrix(i,j) = norm(V-V2);
    end
end     

end

