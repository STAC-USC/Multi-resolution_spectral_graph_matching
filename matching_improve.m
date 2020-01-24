function [ zmatrix ] = matching_improve( G_1, G_2, zmatrix, numbSecVert, beta)
%Find three anchor point pairs to improve matching

% we want to utilize the info of the distance from 3 anchor points to each vertices in the graph
d_1 = distances(G_1);
d_2 = distances(G_2);
zmatrixanchor = zmatrix;
%     anchorInfo = zeros(numbSecVert,2);
% this for loop is useless
for i=1:numbSecVert
    [~,II] = min(zmatrixanchor(:));
    [II_row, II_col] = ind2sub(size(zmatrixanchor),II);
    %translate their similarity by using three anchor point, which are most similar
    %'zmatrixanchor' is useless
    zmatrixanchor = zmatrixanchor - beta*(repmat(d_1(:,II_row),1,size(zmatrixanchor,2))+repmat(d_2(II_col,:),size(zmatrixanchor,1),1));
    zmatrixanchor(:,II_col) = 255*ones(size(zmatrixanchor,1),1); % typo
    zmatrixanchor(II_row,:) = 255*ones(1,size(zmatrixanchor,2)); % typo
%         anchorInfo(i,:) = [II_row, II_col];
%         
%         II_row = anchorInfo(i,1);
%         II_col = anchorInfo(i,2);
    % 'resta' is the distance residue between a pair of anchor points
    resta = beta*(repmat(d_1(:,II_row),1,size(zmatrix,2))-repmat(d_2(II_col,:),size(zmatrix,1),1));
    % the more residue means the corrsponding vertex index is more unreliable. 
    zmatrix = zmatrix + abs(resta);
end

end

