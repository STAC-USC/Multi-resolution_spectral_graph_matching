function [ matrix1histmod,matrix2histmod ] = pair_eigenvector(K, matrix1red, matrix2red, helpyn, valormin )
%% find the k most similar eigenvector in this function 
    % deal with sign ambiguity
    [ Amatrix,matsign ] = sign_ambiguity( K,matrix1red, matrix2red );
    
    % pair up eigenvectors
    unomas = 1;
    % 'matrix1histmod' and 'matrix2histmod' store the K pairs eigenvectors based on their similarity
    matrix1histmod = [];
    matrix2histmod = [];
    % why use k here not k'
    for i=1:K
        [valor,I] = min(Amatrix(:));
        if valor < valormin
            % find the minimal value's location in 'Amatrix'
            [I_row, I_col] = ind2sub(size(Amatrix),I);
            if helpyn == 0
                help = 1;
            else
                help = Amatrix(I_row,I_col);
            end
            % delete the value by filling a huge value in the matrix once we have paired them up
            Amatrix(:,I_col) = 255*ones(K,1);
            Amatrix(I_row,:) = 255*ones(1,K);
            
            % Order the eigenvector and only keep the ones similar
            % 'I_row' and 'I_col' are the indics of eigenvectors we select from graph 'G1' and 'G2' respectively
            matrix1histmod(:,unomas) = matrix1red(:,I_row).*(1/help);
            matrix2histmod(:,unomas) = (matrix2red(:,I_col).*matsign(I_row,I_col)).*(1/help);
            unomas = unomas +1;
        end
    end

end

