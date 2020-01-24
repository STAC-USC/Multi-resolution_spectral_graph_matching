function [Gs1, Gs2, matrix1, matrix2]=kron_gsp_graph_multiresolution(G1,G2,num_levels,param)
% GSP_GRAPH_MULTIRESOLUTION  Compute a multiresolution of graphs
% Modified version of gsp_graph_multiresolution.m 

if nargin < 4
    param = struct;
end
% Sertting default value for every parameters
if ~isfield(param,'sparsify'), param.sparsify = 1; end;
if ~isfield(param,'compute_full_eigen'), param.compute_full_eigen = 0; end;
if ~isfield(param,'sparsify_epsilon')
    param.sparsify_epsilon1 = min(10/sqrt(G1.N),.3);
    param.sparsify_epsilon2 = min(10/sqrt(G2.N),.3); 
end;

if ~isfield(param,'reduction_method')
    param.reduction_method='kron'; 
end

if ~isfield(param,'downsampling_method')
    param.downsampling_method='largest_eigenvector';
end

% if there is no eigenvector, we do a computation
if param.compute_full_eigen
    if (~isfield(G1,'U') || ~isfield(G1,'e') )
        G1=gsp_compute_fourier_basis(G1);
    end
    if (~isfield(G2,'U') || ~isfield(G2,'e') )
        G2=gsp_compute_fourier_basis(G2);
    end
else
    if ~isfield(G1,'lmax')
        G1=gsp_estimate_lmax(G1);
    end
    if ~isfield(G2,'lmax')
        G2=gsp_estimate_lmax(G2);
    end
end

%set up cell for multiresolutions of graphs
Gs1=cell(num_levels+1,1);
Gs2=cell(num_levels+1,1);
Gs1{1}=G1;
Gs1{1}.mr.idx=(1:Gs1{1}.N)';
Gs1{1}.mr.orig_idx=Gs1{1}.mr.idx;
Gs2{1}=G2;
Gs2{1}.mr.idx=(1:Gs2{1}.N)';
Gs2{1}.mr.orig_idx=Gs2{1}.mr.idx;

suma = 2;
for lev=1:num_levels
    Gs1{lev+1}.directed=0;
    Gs2{lev+1}.directed=0;
    
    % Graph downsampling: get indices to keep for the new lower resolution graph
    switch param.downsampling_method
        case 'degree'
            % sum adjaency matrix on second dimension
            degree1 = sum(Gs1{lev}.W,2); 
            % select half of the most densly connected nodes
            [~,indexs1]=sort(full(degree1),'descend');
            % 'keep_inds1' stores the indices that are selected
            keep_inds1 = indexs1(1:round(length(indexs1)/2)); 
            matrix1(1:length(keep_inds1),lev)=keep_inds1;
            
            degree2 = sum(Gs2{lev}.W,2); 
            [~,indexs2]=sort(full(degree2),'descend'); 
            keep_inds2 = indexs2(1:round(length(indexs2)/2)); 
            matrix2(1:length(keep_inds2),lev)=keep_inds2;
        
        case 'degree_gap'
            degree1 = sum(Gs1{lev}.W,2); 
            [ordering1,indexs1]=sort(full(degree1),'descend');                  
            orderingp1 = circshift(ordering1,-1);
            % residue after shifting
            c = ordering1(1:length(ordering1)-1)-orderingp1(1:length(ordering1)-1);
            % find interval index of 0.4 and 0.6
            first = round((length(c)/2)-(0.1*length(c)));
            endf = round((length(c)/2)+(0.1*length(c)));
            % find the index of maximum gap, 'indexc1'
            [~,indexc1] = max(c((first):(endf)));
            indexc1 = indexc1 +(first-1);
            keep_inds1 = indexs1(1:indexc1); 
            matrix1(1:length(keep_inds1),lev)=keep_inds1;   
                        
            degree2 = sum(Gs2{lev}.W,2); 
            [~,indexs2]=sort(full(degree2),'descend');                  
            indexc2 = indexc1;
            keep_inds2 = indexs2(1:indexc2); 
            matrix2(1:length(keep_inds2),lev)=keep_inds2;
            
        case 'degreeOnlyfirst'
            degree1 = sum(Gs1{1}.W,2); 
            [~,indexs1]=sort(full(degree1),'descend'); 
            keep_inds1 = indexs1(1:round(length(indexs1)/suma));             
            [~, indexa] = ismember(Gs1{lev}.coords,Gs1{1}.coords(keep_inds1,:),'rows');
            keep_inds1 = find(indexa > 0);            
            matrix1(1:length(keep_inds1),lev)=keep_inds1;
            
            degree2 = sum(Gs2{1}.W,2); 
            [~,indexs2]=sort(full(degree2),'descend'); 
            keep_inds2 = indexs2(1:round(length(indexs2)/suma));             
            [~, indexb] = ismember(Gs2{lev}.coords,Gs2{1}.coords(keep_inds2,:),'rows');
            keep_inds2 = find(indexb > 0);            
            matrix2(1:length(keep_inds2),lev)=keep_inds2;
            suma = suma*2;
            
       case 'freqbest'
           % 'e' contains eigenvalue
            number = find(Gs1{lev}.e>0.001);
            % select the eigenvector corresponding to the non-zero smallest eigenvalue 
            t1 = Gs1{lev}.U(:,number(1));
            t2 = Gs2{lev}.U(:,number(1));
            t3 = -Gs2{lev}.U(:,number(1));

            [t1,indexs1] = sort(t1);
            t2 = sort(t2);
            t3 = sort(t3);
            
            mit = 0; %(cotsup+cotinf)/2;
            
            % Although it may seem it is the opposite, it is not due to
            % the matlab method sort()
            if sum(abs(t1-t2)) > sum(abs(t1-t3))
                keep_inds2 = find(Gs2{lev}.U(:,number(1))>mit);
            else
                keep_inds2 = find((-Gs2{lev}.U(:,number(1)))>mit);
            end
                 
            keep_inds1 = indexs1(1:length(keep_inds2));
                        
            matrix1(1:length(keep_inds1),lev)=keep_inds1;
            matrix2(1:length(keep_inds2),lev)=keep_inds2;
        
        case 'largest_eigenvector'
            % Based on the polarity of the largest eigenvector
            if isfield(Gs1{lev},'U')
                largest_eigenvector1 = Gs1{lev}.U(:,Gs1{lev}.N);
            else
                [largest_eigenvector1,~]=eigs(Gs1{lev}.L,1); 
            end
            largest_eigenvector1=largest_eigenvector1*sign(largest_eigenvector1(1));
            nonnegative_logicals1=(largest_eigenvector1 >= 0);
            if sum(nonnegative_logicals1) == 0
                error('Too many pyramid levels. Try fewer.');
            end
            
            if lev ==1
                matrix1 = zeros(length(find(nonnegative_logicals1))+1,num_levels); 
            end 
                
            keep_inds1=find(nonnegative_logicals1);
            matrix1(1:length(keep_inds1),lev)=keep_inds1;
            
            if isfield(Gs2{lev},'U')
                largest_eigenvector2 = Gs2{lev}.U(:,Gs2{lev}.N);
            else
                [largest_eigenvector2,~]=eigs(Gs2{lev}.L,1); 
            end
            largest_eigenvector2=largest_eigenvector2*sign(largest_eigenvector2(1));
            nonnegative_logicals2=(largest_eigenvector2 >= 0);
            if sum(nonnegative_logicals2) == 0
                error('Too many pyramid levels. Try fewer.');
            end
            
            if lev ==1
                matrix2 = zeros(length(find(nonnegative_logicals2))+1,num_levels); 
            end 
                
            keep_inds2=find(nonnegative_logicals2);
            keep_inds2 = keep_inds2(1:length(keep_inds1));
            matrix2(1:length(keep_inds2),lev)=keep_inds2;
            
        otherwise
            disp('This graph downsampling method does not exist');
    end
    
   
    % Graph reduction: rewire the new lower resolution graph to form weighted adjacency and Laplacian matrices
    switch param.reduction_method
        case 'kron'
            % Kron reduction
            Gs1{lev+1}.L=gsp_kron_reduce(Gs1{lev}.L,keep_inds1);
            Gs2{lev+1}.L=gsp_kron_reduce(Gs2{lev}.L,keep_inds2);        
        otherwise
            error('Unknown graph reduction method');
    end
 
    % Create the new graph from the reduced weighted adjacency matrix 
    % N is the total number of nodes; L is laplacian matrix; A is sign
    % matrix; W is adjacency matrix;
    Gs1{lev+1}.N=size(Gs1{lev+1}.L,1);
    Gs1{lev+1}.W=diag(diag(Gs1{lev+1}.L))-Gs1{lev+1}.L;
    Gs1{lev+1}.A=sign(Gs1{lev+1}.W);
    Gs1{lev+1}=gsp_copy_graph_attributes(Gs1{lev},1,Gs1{lev+1});
    
    Gs2{lev+1}.N=size(Gs2{lev+1}.L,1);
    Gs2{lev+1}.W=diag(diag(Gs2{lev+1}.L))-Gs2{lev+1}.L;
    Gs2{lev+1}.A=sign(Gs2{lev+1}.W);
    Gs2{lev+1}=gsp_copy_graph_attributes(Gs2{lev},1,Gs2{lev+1});
    
    % Spectral sparsification
    % This section can be deleted since graph sparsification is not
    % performed
    if param.sparsify
        if Gs1{lev+1}.N>2
            sparsify_epsilon=max(param.sparsify_epsilon1,2/sqrt(Gs1{lev+1}.N));
            Gs1{lev+1} = gsp_graph_sparsify(Gs1{lev+1},sparsify_epsilon);
        end
        if Gs2{lev+1}.N>2
            sparsify_epsilon=max(param.sparsify_epsilon2,2/sqrt(Gs2{lev+1}.N));
            Gs2{lev+1} = gsp_graph_sparsify(Gs2{lev+1},sparsify_epsilon);
        end
    else
           Gs1{lev+1}.Ne = nnz(Gs1{lev+1}.W)/2;
           Gs2{lev+1}.Ne = nnz(Gs2{lev+1}.W)/2;
    end
    
    % Copy the coordinates of the subsampled vertices
    Gs1{lev+1}.coords = Gs1{lev}.coords(keep_inds1,:);
    Gs1{lev+1}.type='from multiresolution';
    
    Gs2{lev+1}.coords = Gs2{lev}.coords(keep_inds2,:);
    Gs2{lev+1}.type='from multiresolution';
    
    % Update indexing
    Gs1{lev+1}.mr.idx=keep_inds1;
    Gs1{lev+1}.mr.orig_idx=Gs1{lev}.mr.orig_idx(keep_inds1);
    
    Gs2{lev+1}.mr.idx=keep_inds2;
    Gs2{lev+1}.mr.orig_idx=Gs2{lev}.mr.orig_idx(keep_inds2);
    
    % Compute full eigendecomposition of new graph, if desired
    if param.compute_full_eigen
        Gs1{lev+1}=gsp_compute_fourier_basis(Gs1{lev+1});
        Gs2{lev+1}=gsp_compute_fourier_basis(Gs2{lev+1});
    else
        Gs1{lev+1}=gsp_estimate_lmax(Gs1{lev+1});
        Gs2{lev+1}=gsp_estimate_lmax(Gs2{lev+1});
    end
     
end
end



  


