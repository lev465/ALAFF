% format long
% 
% rng;
% 
% T = rand(1, 1);
% eigen = eig(T);
% T = Spectral_Decomposition_Lambda(T);
% %assert (all(eigen) == all(eig(T)));
% 
% T = rand (2, 2);
% eigen = eig (T);
% T = Spectral_Decomposition_Lambda(T);
% %assert (all (eigen) == all (eig(T)));
% 
% for m = 3:20
%     % Create a random m x m tridiagonal matrix.  Even though we will only
%     % update the lower triangular part, we create the symmetric matrix.
%     m
%     T = rand( m, m );
% 
%     % Extract upper triangle plus first subdiagonal.
%     T = triu( T, -1 );
% 
%     % Set strictly upper triangular part to transpose of strictly lower
%     % triangular part
%     T = tril( T ) + tril( T,-1 )';
% 
%     % Make a copy of T
%     T1 = T;
% 
%     eigen = eig(T1);
% 
%     T1 = Spectral_Decomposition_Lambda (T1)
% 
%     [eigen eig(T1)];
% 
%     assert (all (abs(eig(T1)) - abs(eigen) <=  1e-10))
% 
% end

function T = Spectral_Decomposition_Lambda( T )
    assert (size(T,1) == size(T,2), 'matrix is not square');

    if (size(T,1) == 1)
        return;
    end

    m = size(T,1);
    T1 = T;

    while m > 2
        while abs( T1(m, m - 1) ) > eps(1)* ( sqrt( abs (T1 (m-1, m-1) )+ abs (T1 (m, m) ) ) ) 
            T1 = Francis_Step (T1);
        end

        T1 = tril( T1 ) + tril( T1, -1 )';

        T(m, m) = T1(m, m);

        m = m - 1;

        T1 = T1(1:m, 1:m);
    end

    T1 = eig(T1);

    T (1:2, 1:2) = [T1(1, 1) 0; 0 T1(2,1)];

    T = diag(diag(T));
end