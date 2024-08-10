format long

% matrix size
    for i=1:1000
    m = 5;
    
    % Generate a random matrix
    A = rand( m, m );
    
    % Create vector in which to store the scalars tau from the Householder
    % transformations
    t = rand( m, 1 );
    
    % Create vector in which to store the scalars rho from the Householder
    % transformations
    r = rand( m, 1 );
    
    % Compute the reduction to tridiagonal form
    [ B, t, r ] = BiRed( A, t, r );
    
    % Quick check if it was probably done correctly: Check the singular values of
    % the bidiagonal matrix (extracted from B) with those of the original matrix.
    Bi = BiFromB( B )
    disp('Singular values of B')
    svd( Bi )
    
    disp('Singular values of A')
    svd( A )

    % Bi = [   -1.5990    2.1082         0         0         0;
    %      0   -0.8816    0.5861         0         0;
    %      0         0    0.4670    0.5783         0;
    %      0         0         0   -0.3608    0.1884;
    %      0         0         0         0    0.0757]


    [S, U, V] = SVD_BiDiag_ImpShift(Bi);

    U;
    diag(S)
    V;

    [U_bi, S_bi, V_bi] = svd(Bi);
    
    format short
    U./U_bi;
    diag(S)./S_bi
    assert(all(S - diag(S_bi) < 1e-8))


    assert( all (svd (Bi) - svd (A) < 1e-5), 'different singular values');
    end
