% Implementation of Umeyama Spectral Graph Matching Algorithm
%% Create the two Graphs
%[G1, G2] = kron_createG1G2(1);

%% Obtain U
G1=gsp_compute_fourier_basis(G1);
G2=gsp_compute_fourier_basis(G2);

%% Drawing
param.climits = [0 1];

matrix1red = G1.U;
matrix2red = G2.U;

%% Sign Correction Methods
matrix1red = abs(matrix1red);
matrix2red = abs(matrix2red);

%% Compute Zmatrix
zmatrix = matrix1red*matrix2red';

    
%% Decide the vector auxi based on the maximums of zmatrix
auxi = zeros(1,size(matrix1red,1));
Pmat = zeros(size(matrix1red,1),size(matrix1red,1));
for i=1:min(size(zmatrix))
    [~,I] = max(zmatrix(:));
    [I_row, I_col] = ind2sub(size(zmatrix),I);
    Pmat(I_row,I_col) = 1;
    auxi(I_row) = I_col;
    zmatrix(I_row,:) = zeros(1,size(matrix2red,1));
    zmatrix(:,I_col) = zeros(size(matrix1red,1),1);
end

%% Compute the error generated in the graph matching
error = norm((Pmat*(G2.W)*Pmat'-(G1.W)),'fro');   

%% Decide f1 and f2 for the drawing
mat = G2.coords;
[~,idx] = sort(mat(:,1));
f2(idx) = linspace(0,1,size(G2.U,1));
f1(auxi~= 0)=f2(auxi(auxi~= 0));

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



