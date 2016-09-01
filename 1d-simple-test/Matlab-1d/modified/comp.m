function [l2_U,l2_R]= comp( Nx )

close all;

% Script to solve 1-D compressible stratified waves
% with a Hamiltonian finite volume scheme.

% We assume that the equations have been scaled such that 
% g = c^2_0 = 1 


% Boundaries at x = a & x = L are taken to be solid walls,
% so no normal flow through the boundary. A similar condition is taken 
% for the test function, to preserve skew symmetry, which leads to a
% condition for the density at the boundary.


% Currently attempts to save animation as movies fail.
% Uncomment pause code in the timestep loop to view wave evolution.

% Need to add clause for Nx = 2 , 3 in matrix loops

%% Mesh Constants

a = 0; % Mesh starting point
L = 1; % Mesh end point
%Nx = 64; % Number of elements
dx = (L-a)/Nx; % Element size
dt = dx; % Timestep discretisation

theta0 = 0.33; % Flux constant, 0 < theta < 1
periods = 10; % End time for simulation, as period = 1 non dim time
N_sq = 3;
m = 2 * pi;
sigma = sqrt(0.25 * N_sq ^ 2  + m ^ 2);

t = 0; % Starting time

%% Storage
% Vector for solution
% U = [U,R]^T
U = zeros(2*Nx,1);

R_0_centre = zeros(Nx,1);
R_0_face = zeros(Nx+1,1);

theta = [ 0.5; 0.5; ones(Nx-3,1)*theta0; 0.5; 0.5 ]; % variable theta 

dR_0_centre =  zeros(Nx,1);

% Discrete Energy
H = zeros(1+periods/dt,1);

% Storage for flux matrices
DIV = zeros(2*Nx, 2*Nx);

%% initial data

% cell and face grid
X = linspace(a+dx/2, L-dx/2, Nx)'; % cells   %%% WB: X = linspace(a+dx, L-dx, Nx); 
x_node = linspace(a, L, Nx+1)'; % faces

% background density
R_0_centre = exp( -3.0*X );
dR_0_centre = -3.0*exp( -3.0*X );

R_0_face = exp( -3.0*x_node );

% R_0_centre = 0.5*( R_0_face(1:Nx) + R_0_face(2:Nx+1) );
% dR_0_centre = ( R_0_face(2:Nx+1) - R_0_face(1:Nx) )/dx;

% initial conditions
U(1:Nx) = sin(sigma*0.125)*exp(-0.5*N_sq*X).*sin(m*X); % Velocity
U(Nx+1:2*Nx) = cos(sigma*0.125)*exp(-0.5*N_sq*X).*( 0.5*N_sq*sin(m*X) + m*cos(m*X) )/sigma; % Density

% Initial Energy
H0 = 0.5*dx*sum( ( U(1:Nx).*U(1:Nx) + U(Nx+1:2*Nx).*U(Nx+1:2*Nx) )./R_0_centre );
H(1) = H0;

%% Build flux matrices
% First we construct the system d/dt(U) = DIV(U),
% DIV(U) will be the discrete divergence .
% In this example DIV is the sum of the numerical line fluxes for
% each element.

% Internal cell contributions - this will form a tridiagonal system
for i = 2:Nx-1
    % Flux for velocity U, which depends on R
    DIV(i, Nx + i -1) =  theta(i) / ( dx * R_0_centre(i -1 )) * R_0_face( i ) ;
    DIV(i, Nx + i ) =  -theta(i) / ( dx * R_0_centre(i ))* R_0_face( i ) ...
                   + (1- theta(i+1))/ ( dx * R_0_centre(i  ))* R_0_face( i+1 ) ...
                   - dR_0_centre(i)/R_0_centre(i);
    DIV(i, Nx + i +1) =  -(1 - theta(i+1))/ ( dx * R_0_centre(i + 1 )) * R_0_face( i + 1 ) ; 
    % Flux for density R, which depends on U
    DIV(Nx + i,  i -1) = (1 - theta(i))/ ( dx * R_0_centre(i -1 )) * R_0_face( i ) ;
    DIV(Nx + i,  i ) =  theta(i)/ ( dx * R_0_centre(i  ))* R_0_face( i ) ...
                  - (1-theta(i+1))/ ( dx * R_0_centre(i ))* R_0_face( i+1 ) ...
                  +  dR_0_centre(i)/R_0_centre(i);
    DIV(Nx + i,  i +1) = - theta(i+1)/ ( dx * R_0_centre(i + 1 )) * R_0_face( i + 1 ) ; 
end    

% Solid wall boundaries
% Left wall
DIV(1,Nx+1) =  (1-theta(2))/ ( dx * R_0_centre(1 )) * R_0_face( 2 ) -  dR_0_centre(1)/R_0_centre(1);
DIV(1,Nx+2) = - (1-theta(2))/ ( dx * R_0_centre(2 )) * R_0_face( 2 );
DIV(Nx+1,1) = -(1-theta(2))/ ( dx * R_0_centre(1 )) * R_0_face( 2 ) +  dR_0_centre(1)/R_0_centre(1);
DIV(Nx+1,2) = -theta(2)/ ( dx * R_0_centre(2 )) * R_0_face( 2 );
% Right wall
DIV(Nx,2*Nx-1) = theta(Nx)/ ( dx * R_0_centre(Nx -1  )) * R_0_face( Nx );
DIV(Nx,2*Nx) = -theta(Nx)/ ( dx * R_0_centre(Nx )) * R_0_face( Nx ) -  dR_0_centre(Nx)/R_0_centre(Nx);
DIV(2*Nx,Nx-1) = (1-theta(Nx))/ ( dx * R_0_centre(Nx -1  )) * R_0_face( Nx );
DIV(2*Nx,Nx) = theta(Nx)/ ( dx * R_0_centre(Nx )) * R_0_face( Nx ) +  dR_0_centre(Nx)/R_0_centre(Nx);

% Create timestepper matrices
% Using implicit mid-point timestepper creates matrices acting on 
% U^(n+1) and U^n.
% The system we will solve will be 
% P U^(n+1) = Q U^n
% P = I - 0.5*dt*DIV
% Q = I + 0.5*dt*DIV

P = eye(2*Nx) - (dt)*0.5*DIV;
Q = eye(2*Nx) + (dt)*0.5*DIV;

% Inverse of P\Q is not time dependant so we will build it here
%Inverse = P\Q;

[low,up,piv] = lu(P);

%F(ceil(periods/dt)) = struct('cdata',[],'colormap',[]);

%% Timestep Loop
count_energy = 1;
figure

tdone=0; num=0;

while t < periods
    
    t =t+dt;
    
    tdone = tdone+dt;
    if( tdone > 1 )
        num = num+1;
        fprintf('Period %d\n',num);
        tdone=0;
    end
    
    % Calculate  U^(n+1) from  U^n and advance in time
    U=solvePQ(low,up,piv,P,Q*U);  %% WB: U = Inverse*U;
    %U=solvePQ(P,Q*U);
   
    % Calculate Energy
    count_energy = count_energy+1;
    H(count_energy,1) = 0.5*dx*sum( (U(1:Nx).*U(1:Nx) + U(Nx+1:2*Nx).*U(Nx+1:2*Nx))./R_0_centre );
    
    % Plot velocity 
    plot(X,U(1:Nx),'-b',X,U(Nx+1:2*Nx),'-r')
    title(['Solution after ', num2str(dt*count_energy), ' periods - Velocity u'])
    axis([  a L -1 1])
    xlabel('x')
    ylabel('u, rho')
    %legend('u','rho')
    %pause(0.01)
    drawnow
    hold off
    
end

%F(count_energy)=getframe(gcf);
%movie2avi(F, 'velocity.avi')

% Plot Error in Energy
figure
hold all
set(gca,'fontsize', 12)
plot(linspace(0,t,length(H)), (H-H0)/H0); %%% WB plot(linspace(0,periods,periods/dt), (H-H0)/H0)
title('Relative Error in Energy')
xlabel('Number of Periods')
hold off

% Exact Solution at time = t_end
U_exact(1:Nx,1) = exp(-0.5 * N_sq *X) .* sin(m*X) * sin(sigma*(t+0.125)); % Velocity
U_exact(Nx+1:2*Nx,1) = exp(-0.5 * N_sq *X) .* cos(sigma*(t+0.125)) .* ( 0.5*N_sq*sin(m*X) + m*cos(m*X) )/sigma; % Density

error_U = sum( (U(1:Nx)-U_exact(1:Nx)).^2 );
error_R = sum( (U(Nx+1:2*Nx)-U_exact(Nx+1:2*Nx)).^2 );

% 
% l_infty_U = max(error_U);
% l_infty_R = max(error_R);
l2_U = sqrt(sum(error_U)/Nx);
l2_R = sqrt(sum(error_R)/Nx);

end


%% 
% extended precision linear solve for Ax=b
% A=(L,U,P) factorisation

function x = solvePQ( L,U,P, A,b)

x = P*( U\( L\b ) );          % solve A x = b

r = double(b-A*sym(x,'f'));   % compute residual r = b - A x

x = x + P*( U\( L\r ) );      % solve A e = r and correct x = x + e

end
