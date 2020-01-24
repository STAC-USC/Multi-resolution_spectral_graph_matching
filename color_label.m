function [ f1,f2 ] = color_label( Gs1,matrix1red,auxi,ii )
% 'coords' contains 3D location of all vertices. We transform them and plot them in 2D space
mat = Gs1{numel(Gs1)-ii+1}.coords;
% get the order based on the first coordinate.
[~,idx] = sort(mat(:,1));
% 'f1' and 'f2' are color label
f1(idx) = linspace(0,1,size(matrix1red,1));
f2(auxi~= 0)=f1(auxi(auxi~= 0));
% make the row vector become the column one
if size(f1,1) == 1
    f1 = f1';
end
if size(f2,1) == 1
    f2 = f2';
end
end

