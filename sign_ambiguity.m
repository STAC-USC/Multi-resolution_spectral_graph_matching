function [ Amatrix,matsign ] = sign_ambiguity( K,matrix1red, matrix2red )
%% deal with sign ambiguity
% 'Amatrix' is H matrix descirbed in the paper, which measures the similarity between all pairs of eigenvectors
% 'matsign' store the sign providing the best match for each pair
Amatrix = 255*ones(K,K);
matsign = ones(K,K);
for i=1:K
    for j=1:K
        % select the first K eigenvectors from each graph to compare their similarity
        l1 = sort(matrix1red(:,i));
        l2 = sort(matrix2red(:,j));
        l3 = sort(-matrix2red(:,j));
        % similarity function defined in the paper    
        Amatrix(i,j) = min(sum(abs(l1-l2)),sum(abs(l1-l3)));
        % deal with sign ambiguity by checking which is the most similar
        if sum(abs(l1-l2)) > sum(abs(l1-l3))
            matsign(i,j) = -1;
        end
    end
end

end

