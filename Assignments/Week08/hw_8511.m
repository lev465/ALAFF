close all;

% Converts a sparse matrix into compressed row storage.
%
% Inputs:     [N] m x n Sparse matrix to be deconstructed
%
% Outputs:    [nzA] vector that stores the nonzero elements of matrix N, 
%                     size nnzeroes
%             [ir] vector that stores the index in array nzA where the first 
%             element of the each row is stored, size n+1
% 
%             [ic] vector that holds the column indices of the corresponding 
%             elements in array, size nnzeroes

function [ nzA, ir, ic ] = Create_Poisson_problem_nzA( N )
    nnzeroes = nnz(N);

    nzA = [];
    ir = zeros(size(N, 1) + 1, 1);
    ic = [];


    for i=1:size(N, 1)
        curr_elements = 0;

        for j=1:size(N,2)
            if(N (i,j) ~= 0)
                nzA (end + 1, 1) = N (i,j);
                ic (end + 1, 1) = j;

                if curr_elements == 0
                    ir (i) = size(nzA, 1);
                end

               curr_elements = curr_elements + 1;
            end
        end
    end

    ir (size(N, 1) + 1) = nnzeroes + 1;

    for i=(size(ir, 1) - 1):-1:1
        if (ir (i) == 0)
            ir (i) = ir (i+1);
        end
    end

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


% might be wise to write some edge test cases here
% zero matrix
A = [ 0 0 0 0 0 0; 0 0 0 0 0 0 ; 0 0 0 0 0 0];
x = [1; 2; 3; 4; 5; 6];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% identity
A = [ 1 0 0 0 0 0; 0 1 0 0 0 0 ; 0 0 1 0 0 0];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% normal matrices
A = [0 1 2 3 0 4; 5 1 2 3 5 4; 0 1 2 3 9 4 ];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

A = [0 0 0 0 0 0; 0 1 0 1 0 1; 0 0 0 0 0 0];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% zero vector
A = [0 1 2 3 0 4; 5 1 2 3 5 4; 0 1 2 3 9 4 ];
x = [0; 0; 0; 0; 0; 0];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% zero matrix, zero vector
A = [ 0 0 0 0 0 0; 0 0 0 0 0 0 ; 0 0 0 0 0 0];
x = [0; 0; 0; 0; 0; 0];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% A is a 1-D matrix
A = [0 1 2 3 0 4];
x = [1; 2; 3; 4; 5; 6];
[test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (A);
y = SparseMvMult(test_nzA, test_ir, test_ic, x);
assert (all(y == A*x));

% Randomized testing
for i=1:200
    test_A = randi(100, 8, 20);
    x_test = randi(100, 20, 1);
    [test_nzA, test_ir, test_ic] = Create_Poisson_problem_nzA (test_A);
    y = SparseMvMult(test_nzA, test_ir, test_ic, x_test);
    y_h = test_A*x_test;
    assert(all(y_h == y));
end
