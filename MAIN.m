% ENTREGABLE
% This function performs a graph matching algorithm explained in:
% "Spectral Graph Matching through a Multi-Resolution Approach" by
% Victor Gonzalez Navarro and Antonio Ortega

% To use this method you need to download and install the following library:
% https://epfl-lts2.github.io/gspbox-html/ to plot the graphs

% ------------------------------------------------------------------------------------------------
% Example of input variables for TESTgraph1.mat: 
% [ngraph=1, Nlevel=3, downmethod=3, K=3, valormin=3, numbSecVert=3, beta=0, helpyn=0, pointreg=0]
% ------------------------------------------------------------------------------------------------
% Example of input variables for TESTgraph2.mat: 
% [ngraph=2, Nlevel=3, downmethod=2, K=15, valormin=3, numbSecVert=3, beta=0, helpyn=0, pointreg=0]
% ------------------------------------------------------------------------------------------------
% Example of input variables for TESTgraph3.mat: 
% [ngraph=3, Nlevel=4, downmethod=2, K=15, valormin=6, numbSecVert=3, beta=0, helpyn=0, pointreg=0]
% [ngraph=3, Nlevel=2, downmethod=2, K=10, valormin=3, numbSecVert=3, beta=0, helpyn=0, pointreg=0]
% ------------------------------------------------------------------------------------------------
% Example of input variables for TESTgraph4.mat: 
% [ngraph=4, Nlevel=1, downmethod=2, K=15, valormin=6, numbSecVert=3, beta=0, helpyn=0, pointreg=0]
% ------------------------------------------------------------------------------------------------

%% General variables
Nlevel=3; 
% graph downsampling method
downmethod=2;
% the value of K to select the first K eigenfunctions
K=3;
% the number of secure vertices
valormin=3;
% the number of secure vertices
numbSecVert=3;
% the value of beta (importance due to secure vertices)
beta=0;
% whether to use the modified euclidean distance
helpyn=0;
% whether to perform a point registragion
pointreg=0;

param.sparsify = 0;
param.compute_full_eigen = 1;
if downmethod == 1
    param.downsampling_method = 'degree_gap'; %using the gap method
elseif downmethod == 2
    param.downsampling_method = 'freqbest'; %best for pointclouds
elseif downmethod == 3
    param.downsampling_method = 'degree'; %best for random graphs
end

%% apply graph downsampling and kron reduction to obtain the multiresolution 
% 'Gs1' and 'Gs2' are new graphs of 'G1' and 'G2' with 'Nlevel' resolution
[ Gs1, Gs2, keep_inds1, keep_inds2 ] = kron_gsp_graph_multiresolution(G1, G2, Nlevel, param);

%% Drawing
param.climits = [0 1];
param.show_edges = 0;
param.colorbar = 0;
figure;

%% graph matching based on multiresolution
for ii = 1:numel(Gs1)
    K = K+1;
    % 'matrix1red' and 'matrix2red' contain all eigenvectors of G1 and G2 repectively
    matrix1red = Gs1{numel(Gs1)-ii+1}.U;
    matrix2red = Gs2{numel(Gs2)-ii+1}.U;
    
    % Find matching matrix in step4
    [ zmatrix,matrix1histmod,matrix2histmod ] = matching_matrix(K, matrix1red, matrix2red, helpyn, valormin );

    % Apply point registration algorithm 
    if pointreg == 1
        [RotMat,TransVec,matrix2histmod]=icp_pr(matrix1histmod',matrix2histmod');
        matrix2histmod = matrix2histmod'; 
    end
    
    % Find the three anchor point pairs to improve matching
    G_1 = graph(Gs1{numel(Gs1)-ii+1}.W);
    G_2 = graph(Gs2{numel(Gs2)-ii+1}.W);
    [ zmatrix ] = matching_improve(G_1, G_2, zmatrix, numbSecVert, beta);
    
    % find permutation matrix and indices pairs based on matching matrix, zmatrix
    [ Pmat, auxi ] = permutation_matrix(zmatrix, matrix1red, matrix2red);
    
    % Compute the error, J(P) in paper, generated in the graph matching of every resolution
    error = norm((Pmat*(Gs2{numel(Gs1)-ii+1}.W)*Pmat'-(Gs1{numel(Gs1)-ii+1}.W)),'fro');
    
    % Algorithm to compute the final zmatrix based on the zmatrix of every resolution and the error
    if ii > 1
        %-------- Eliminate Pmat and maybe its better----------
        zmatrix = (zmatrix*(error+1))-(Pmat);
        % flip then shift 'zmatrixAnt'
        zmatrixAnt = max(max(zmatrixAnt)) - zmatrixAnt;
        % the way to aggregate the matching matrix of all resolution
        zmatrix(keepinds1,keepinds2) = zmatrix(keepinds1,keepinds2) - zmatrixAnt;
        % pass value to 'zmatrixAnt' for next resolution
        zmatrixAnt = zmatrix;
    else
        % 'zmatrixAnt' store the similarity info of the last resolution
        zmatrixAnt = zmatrix*(error+1); %-(Pmat);
    end
    
    % find the color label. 'f1' and 'f2' are color label for each graph.
    [ f1,f2 ] = color_label( Gs1,matrix1red,auxi,ii );
    
    subplot(2,Nlevel+2,-ii+Nlevel+2)
    zkron_gsp_plot_signal(Gs1{numel(Gs1)-ii+1},f1,param,1)
    title(['\alpha = J(P) = ',num2str(round(error,2))])
    xlabel(['G1 res. level: ', num2str(Nlevel-ii+2)]);

    subplot(2,Nlevel+2,-ii+Nlevel+7)
    zkron_gsp_plot_signal(Gs2{numel(Gs1)-ii+1},f2,param,1)
    xlabel(['G2 level: ', num2str(ii-1)]);
    xlabel(['G2 res. level: ', num2str(Nlevel-ii+2)]);
    
    % Find the indices kept for the next resolution
    if ii<(Nlevel+1)
        keepinds1 = keep_inds1(keep_inds1(:,numel(Gs1)-ii) ~= 0, numel(Gs1)-ii);
        keepinds2 = keep_inds2(keep_inds2(:,numel(Gs2)-ii) ~= 0, numel(Gs2)-ii);               
    end
end

%% Plot and find the final labels
[ Pmat, auxi] = permutation_matrix( zmatrix, matrix1red, matrix2red );
% final error
error = norm((Pmat*(Gs2{numel(Gs1)-numel(Gs1)+1}.W)*Pmat'-(Gs1{numel(Gs1)-numel(Gs1)+1}.W)),'fro');
[ f1,f2 ] = color_label( Gs1,matrix1red,auxi,numel(Gs1) );

figure;
subplot(2,1,1)
zkron_gsp_plot_signal(Gs1{numel(Gs1)-ii+1},f1,param,0)
title(['Final Error = J(P) = ',num2str(round(error,2))])
ax = gca;
ax.FontSize = 15;

subplot(2,1,2)
zkron_gsp_plot_signal(Gs2{numel(Gs1)-ii+1},f2,param,0)


