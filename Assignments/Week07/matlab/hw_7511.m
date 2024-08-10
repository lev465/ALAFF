close all

MAX_ITERATIONS = 10000;
curr_err = 1e-2;
results = zeros(3, 10);
for i = 1:10
    results(1, i) = curr_err;
    [results(2, i), results(3, i)] = GS_and_SOR(MAX_ITERATIONS, curr_err);
    curr_err = curr_err/10;
end
x = string(results(1, :));
b = bar(x, [results(3, :); results(2, :)], 0.6)
legend('SOR','Gauss-Seidel', Location='northwest')
title('Gauss-Seidel vs SOR Iterative Solutions')
xlabel('Maximum Error');
ylabel('Iterations')
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

function [k_GS, k_SOR] = GS_and_SOR(max_iter, max_err)
    
    % Set number of iterations to be performed
    nk = max_iter;
    
    % Set the error for convergence
    max_err = max_err;
    
    % Set parameters alpha and beta
    alpha = 2;
    beta  = 3;
    
    % Set the number of meshpoints so the interior has N x N such points
    N = 50;
    
    % Compute the distance between mesh points, in each direction
    h = 1/(N+1);
    
    % Optimal relaxation parameter - https://userpages.umbc.edu/~gobbert/papers/YangGobbert2007SOR.pdf
    w = 2 / (1 + sin (pi * h));
    
    % We will have arrays that capture the boundary as well as the interior
    % meshpoints.  As a result, we need those arrays to be of size (N+2) x
    % (N+2) so that those indexed 2:N+1, 2:N+1 represent the interior.  
    
    % Compute the x-values at each point i,j, including the boundary
    x = h * [ 0:N+1 ];   % Notice this creates a row vector
    
    % Compute the y-values at each point i,j, including the boundary
    y = h * [ 0:N+1 ];   % Notice this creates a row vector
    
    % Create an array that captures the load at each point i,j
    for i=1:N+2
        for j=1:N+2
            F( i,j ) = ...
                ( alpha^2 + beta^2 ) * pi^2 * sin( alpha * pi * x( i ) ) * sin( beta * pi * y( j ) );
        end
    end
    
    % Set the initial values at the mesh points
    U = zeros( N+2, N+2 );
    
    % Set base error
    err = F - U;
    k_GS = 0;
    
    % GS
    while (sum(abs(err) <= max_err, "all") ~= (N+2)^2)
        U_old = U;
        % update all the interior points
        for i=2:N+1
            for j=2:N+1
                U( i,j ) = ( U( i, j-1 ) + U( i-1, j ) + U( i+1, j ) + U( i, j+1 ) + h^2 * F( i, j ) ) / 4;
            end
        end 
        %mesh( x, y, U );
        %axis( [ 0 1 0 1 -1.5 1.5 ]);
    
        err = U - U_old;
    
        k_GS = k_GS+1;
    
        if k_GS > nk
            break;
        end
        % wait to continue to the next iteration
        %next = input( 'press RETURN to continue' );
    end
    
    k_GS = k_GS-1; % prints last iteration
    
    U_SOR = zeros( N+2, N+2);
    
    % Reset error and counter
    err = F - U_SOR;
    k_SOR = 0;
    
    % SOR
    while (sum(abs(err) <= max_err, "all") ~= (N+2)^2)
        U_old = U_SOR;
        % update all the interior points
        for i=2:N+1
            for j=2:N+1
                %SOR algorithm change
                U_SOR( i,j ) =  ( ( ( U_SOR( i, j-1 ) + U_SOR( i-1, j ) + ...
                    U_SOR( i+1, j ) + U_SOR( i, j+1 ) + h^2 * F( i, j ) ) / 4) + (1 - w) * U_SOR( i,j ) / w) * w;
            end
        end 
        %mesh( x, y, U );
        %axis( [ 0 1 0 1 -1.5 1.5 ]);
        err = U_SOR - U_old;
        k_SOR = k_SOR+1;
    
        if k_SOR > nk
            break;
        end
        % wait to continue to the next iteration
        %next = input( 'press RETURN to continue' );
    end
    k_SOR = k_SOR-1;
end

