% Implementation of David Knossow et al. Spectral Graph Matching Algorithm
%% Create the two Graphs
%[G1, G2] = kron_createG1G2(1);

%% Obtain U
G1=gsp_compute_fourier_basis(G1);
G2=gsp_compute_fourier_basis(G2);

%% Drawing
param.climits = [0 1];
param.show_edges = 0;
param.colorbar = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random Graph
% K = 15;
% valormin = 30;
% Nbins = 21;

%% Bunny
% K = 15;
% valormin = 1000;
% Nbins = 61;

%% Dragon
% K = 15;
% valormin = 1000;
% Nbins = 61;

%% Human
K = 15;
%valormin = 1000;
Nbins = 61;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix1red = G1.U(:,1:K);
matrix2red = G2.U(:,1:K);
matrix1red = normr(matrix1red);
matrix2red = normr(matrix2red);

Amatrix = 255*ones(K,K);
matsign = ones(K,K);
    
%% Compute Amatrix, el kit de la cuestion (PREPROCESSING)
for i=1:K
    for j=1:K
        l1 = hist(matrix1red(:,i),Nbins);
        l2 = hist(matrix2red(:,j),Nbins);
        l3 = hist((-matrix2red(:,j)),Nbins);

        Amatrix(i,j) = min(sum(abs(l1-l2)),sum(abs(l1-l3)));
        if sum(abs(l1-l2)) > sum(abs(l1-l3))
            matsign(i,j) = -1;
        end
    end
end

%% Find matrix1histmod and matrix2histmod
unomas = 1;
matrix1histmod = [];
matrix2histmod = [];
for i=1:K
    [valor,I] = min(Amatrix(:));
    %if valor < valormin
        [I_row, I_col] = ind2sub(size(Amatrix),I);
        Amatrix(:,I_col) = 255*ones(K,1);
        Amatrix(I_row,:) = 255*ones(1,K);

        % Order the eigenvector and only keep the ones similar
        matrix1histmod(:,unomas) = matrix1red(:,I_row);
        matrix2histmod(:,unomas) = (matrix2red(:,I_col).*matsign(I_row,I_col));
        unomas = unomas +1;
    %end
end
unomas

%% Find Zmatrix, using the theory of abs(U)abs(U2) = I or not
zmatrix = 255*ones(size(matrix1histmod,1),size(matrix2histmod,1));
for i=1:size(matrix1histmod,1)
    V = matrix1histmod(i,:);
    for j=1:size(matrix2histmod,1)
        V2 = matrix2histmod(j,:);
        zmatrix(i,j) = norm(V-V2);
    end
end

zmatrixp = zmatrix;

%% Decide the vector auxi based on the minimum of zmatrix
auxi = zeros(1,max(size(zmatrix)));
Pmat = zeros(size(matrix1red,1),size(matrix2red,1));
for i=1:max(size(zmatrix))
    [~,I] = min(zmatrix(:));
    [I_row, I_col] = ind2sub(size(zmatrix),I);
    Pmat(I_row,I_col) = 1;
    auxi(I_col) = I_row;
    zmatrix(:,I_col) = 255*ones(size(zmatrix,1),1);
    zmatrix(I_row,:) = 255*ones(1,size(zmatrix,2));
end

%% Compute the error generated in the graph matching
error = norm((Pmat*(G2.W)*Pmat'-(G1.W)),'fro'); 

%% Decide f1 and f2 for the drawing
mat = G1.coords;
[~,idx] = sort(mat(:,1));
f1(idx) = linspace(0,1,size(matrix1red,1));
f2(auxi~= 0)=f1(auxi(auxi~= 0));

if size(f1,1) == 1
    f1 = f1';
end
if size(f2,1) == 1
    f2 = f2';
end

%% Plot the Graphs
figure;
subplot(2,1,1)
zkron_gsp_plot_signal(G1,f1,param,0)
title(['Error = ',num2str(error)])

subplot(2,1,2)
zkron_gsp_plot_signal(G2,f2,param,0)