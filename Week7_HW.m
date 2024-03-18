close all

% Set number of iterations to be performed
nk = 500;

% Set parameters alpha and beta
alpha = 2;
beta  = 3;

% Set the number of meshpoints so the interior has N x N such points
N = 50;

% Compute the distance between mesh points, in each direction
h = 1/(N+1);
%w = 1.8;
w = 2 / (1 + sin(h*pi)); %Omega


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

u = zeros( N+2, N+2 );
omega = 2 / (1 + sqrt(2 - 2*cos(h*pi)));
%To find the exact solution using SSOR
 
for k = 2 : 1000
        % Natural ordering
        for j = 2 : (N + 1)
            for i = 2 : (N + 1)
                u(i, j) = (1 - omega) * u(i, j) + omega * (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h^2*F(i, j)) / 4;
            end
        end
        
        % Reverse ordering
        for j = (N + 1) : -1 : 2
            for i = (N + 1) : -1 : 2
                u(i, j) = (1 - omega) * u(i, j) + omega * (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) + h^2*F(i, j)) / 4;
            end
        end
end


% Set the initial values at the mesh points.  Use J to indicate
% the array for the Jacobi iteration.  Use GS to indicate the array
% is for the Gauss-Seidel iteration.
J = zeros( N+2, N+2 );
GS = zeros( N+2, N+2 );
SOR  = zeros( N+2, N+2 );
SSOR  = zeros( N+2, N+2 );

error_GS =  [ 0:nk ]*0;
error_SOR =  [ 0:nk ]*0;
error_SSOR =  [ 0:nk ]*0;
error_J =  [ 0:nk ]*0;


% Perform nk iterations
for k = 2:nk
    k           % print current iteration index
    Jold = J;   % we want to use the old values
    
    % update all the interior points (Jacobi iteration)
    for i=2:N+1
        for j=2:N+1
            J( i,j ) = ( Jold( i, j-1 ) + Jold( i-1, j ) + Jold( i+1, j ) + Jold( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end 
    


    % update all the interior points (GS iteration)
    for i=2:N+1
        for j=2:N+1
            GS( i,j ) = ( GS( i, j-1 ) + GS( i-1, j ) + GS( i+1, j ) + GS( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end 

     % update all the interior points (SOR iteration)
    for i=2:N+1
        for j=2:N+1
            SOR( i,j ) = (1-w)*SOR( i, j ) + w*( SOR( i, j-1 ) + SOR( i-1, j ) + SOR( i+1, j ) + SOR( i, j+1 ) + h^2 * F( i, j ) ) / 4;
        end
    end 

    % update all the interior points (SSOR iteration)
        % Forward
        for j = 2 : (N + 1)
            for i = 2 : (N + 1)
                SSOR(i, j) = (1 - omega) * SSOR(i, j) + omega * (SSOR(i - 1, j) + SSOR(i + 1, j) + SSOR(i, j - 1) + SSOR(i, j + 1) + h^2*F(i, j)) / 4;
            end
        end
        
        % Reverse
        for j = (N + 1) : -1 : 2
            for i = (N + 1) : -1 : 2
                SSOR(i, j) = (1 - omega) * SSOR(i, j) + omega * (SSOR(i - 1, j) + SSOR(i + 1, j) + SSOR(i, j - 1) + SSOR(i, j + 1) + h^2*F(i, j)) / 4;
            end
        end
    

    %{
    subplot( 3, 1, 1 );  % plot in top graph
    mesh( x, y, J );
    axis( [ 0 1 0 1 -1.5 1.5 ]);
    
    subplot( 3, 1, 2 );  % plot in bottom graph
    mesh( x, y, GS );
    axis( [ 0 1 0 1 -1.5 1.5 ]);

    subplot( 3, 1, 3);
    mesh( x, y, GS - J );
    axis( [ 0 1 0 1 -0.05 0.05 ])


    subplot( 3, 1, 2);
    mesh( x, y, SOR );
    axis( [ 0 1 0 1 -1.5 1.5 ])
    %}
    if(k>1)
    error_GS(k) = norm(GS -u ) ;
    error_J(k) = norm( J - u) ;
    error_SOR(k) = norm( SOR - u, 'fro') ;   
    error_SSOR(k) = norm(SSOR - u, 'fro') ;
   
    end 

    % wait to continue to the next iteration
   % next = input( 'press RETURN to continue' );
end


%Convergence Plot 
error_GS = log(error_GS);
error_J = log(error_J);
error_SOR = log(error_SOR);

plot(error_GS);
hold on;
plot(error_J);
hold on;
plot(error_SOR);
xlabel('Iterations');
ylabel('Error (log scale)');
legend('GS', 'Jacobi','SOR');
hold off;