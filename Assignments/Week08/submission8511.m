close all;
% ALAFF: Homework 8.5.1.1
% Submission by: Luis Vazquez, LEV465

% Converts a sparse matrix into compressed row storage.
%
% Inputs:     [N] positive integer that denotes size of the 2D Poisson
%                   matrix
%
% Outputs:    [nzA] vector that stores the nonzero elements of matrix N, 
%                     size nnzeroes
%             [ir] vector that stores the index in array nzA where the first 
%             element of the each row is stored, size n+1
% 
%             [ic] vector that holds the column indices of the corresponding 
%             elements in array, size nnzeroes

function [ nzA, ir, ic ] = Create_Poisson_problem_nzA( N )
    D_nnzeroes = N*(1 + 3 * (N-1));
    I_nnzeroes = N * 2 * (N-1);

    nnzeroes = D_nnzeroes + I_nnzeroes;

    nzA = zeros(nnzeroes,1);
    ir = zeros(N^2 + 1, 1);
    ic = zeros(nnzeroes, 1);
    cnt = 1;

    for i=1:N^2
        first_row_filled = false;

        %left-most element
        if (i-N > 0)
            nzA (cnt, 1) = -1;
            ir (i, 1) = cnt;
            ic (cnt, 1) = i-N;
            cnt = cnt + 1;
            first_row_filled = true;
        end

        % add element left diagonal
        if (mod(i-1, N) > 0)
            nzA (cnt, 1) = -1;
            ic (cnt, 1) = i - 1;
            if(~first_row_filled)
                ir (i, 1) = cnt;
                first_row_filled = true;
            end
            cnt = cnt + 1;
        end

        % add diagonal element
        nzA (cnt, 1) = 4;
        ic (cnt, 1) = i;
        if(~first_row_filled)
                ir (i, 1) = cnt;
        end

        cnt = cnt + 1;

        % add right element

        if (mod (i, N) > 0)
            nzA (cnt, 1) = -1;
            ic (cnt, 1) = i + 1;
            cnt = cnt + 1;
        end

        if ( i + N <= N^2)
            
            % add right-most element -I
            nzA (cnt, 1) = -1;
            ic (cnt, 1) = i + N;
            cnt = cnt + 1;

        end

    end

    ir (N^2 + 1) = nnzeroes + 1;


end

% computes y = Ax with the matrix A stored in the sparse format
%
% Inputs:     [nzA] vector that stores the nonzero elements of matrix N, 
%                     size nnzeroes
%             [ir] vector that stores the index in array nzA where the first 
%             element of the each row is stored, size n+1
% 
%             [ic] vector that holds the column indices of the corresponding 
%             elements in array, size nnzeroes
%               
%             [x] vector used for the A*x problem
%
% Outputs:    [y] resulting vector of the A*x problem, where A is the
%               deconstructed matrix

function y = SparseMvMult( nzA, ir, ic, x )
    y = zeros(size(ir, 1) - 1, 1);

    for i = 1:size(y, 1)
        values_in_row = ir(i+1) - ir(i);

        if values_in_row > 0
            temp_row = nzA (ir(i) : ir(i+1) - 1); 
            temp_col = ic (ir(i) : ir(i+1) - 1);
             for j = 1: values_in_row
                y(i) = y(i) + temp_row(j)*x(temp_col(j));
             end
        end
        
    end

end

for i=1:10
    [nzA,ir,ic] = Create_Poisson_problem_nzA(i);
    x = randi(20, i^2, 1);
    y = SparseMvMult(nzA, ir, ic, x);
    %A = makePoissonPdeMatrix(i);
    %assert(all(y) == all(A*x))
end
% 
% function A = makePoissonPdeMatrix(N)
%     A = zeros(N^2);
% 
%     for i=1:N^2
%         %left-most element
%         if (i-N > 0)
%             A(i, i-N) = -1;
%         end
% 
%         % add element left diagonal
%         if (mod(i-1, N) > 0)
%             A(i, i-1) = -1;
%         end
% 
%         % add diagonal element
%         A (i, i) = 4;
% 
%         % add right element
% 
%         if (mod (i, N) > 0)
%             A (i, i+1) = -1;
%         end
% 
%         if ( i + N <= N^2)
% 
%             % add right-most element -I
%             A(i, i+N) = -1;
% 
%         end
% 
%     end
% end