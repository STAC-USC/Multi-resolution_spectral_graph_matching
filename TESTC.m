% Implementation of Terry Caelli and Serhiy Kosinov Spectral Graph Matching Algorithm
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
%valormin = no1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

matrix1red = G1.U(:,1:K);
matrix2red = G2.U(:,1:K);
matrix1red = normr(matrix1red);
matrix2red = normr(matrix2red);
    
%% Compute Amatrix, el kit de la cuestion (PREPROCESSING)
long = length(matrix1red(:,1));
for i=1:K
    l1 = sum(matrix1red(:,i)>0);
    l2 = sum(matrix2red(:,i)>0);
    if (long-l1)>l1
        matrix1red(:,i)=-matrix1red(:,i);
    end
    if (long-l2)>l2
        matrix2red(:,i)=-matrix2red(:,i);
    end
end


%% Find Zmatrix, using the theory of abs(U)abs(U2) = I or not
zmatrix = 255*ones(size(matrix1red,1),size(matrix2red,1));
for i=1:size(matrix1red,1)
    V = matrix1red(i,:);
    for j=1:size(matrix2red,1)
        V2 = matrix2red(j,:);
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