clc; close all; clear all;

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
Nx = 5; % Number of elements
dx = (L-a)/Nx; % Element size
dt = dx; % Timestep discretisation
theta = 0.25; % Flux constant, 0 < theta < 1
periods = 5; % End time for simulation, as period = 1 non dim time
t = 0; % Starting time
N_sq = 3;
m = 2 * pi;
sigma = sqrt(0.25 * N_sq ^ 2  + m ^ 2)


%% Storage
% Vector for solution
% U = [U,R]^T
U = zeros(2*Nx,1);

R_0_centre = zeros(Nx,1);
R_0_face = zeros(Nx+1,1);

dR_0_centre =  zeros(Nx,1);

% Discrete Energy
H = zeros(ceil(periods/dt),1);

% Storage for flux matrices
DIV = zeros(2*Nx, 2*Nx);

%% Initial Condition
X = linspace(a+dx, L-dx, Nx);
x_node = linspace(a, L, Nx+1);

R_0_centre = exp( - 3.0 * X);
R_0_face = exp( - 3.0 * x_node);

dR_0_centre = - 3.0 * exp( - 3.0 * X);

H0 = 0;
% For each element we use the cell centre of the initial condition
for i = 1:Nx
    U(i) = exp(-0.5 * N_sq *X(i)) * sin(m*X(i))*sin(sigma*0.125); % Velocity
    U(Nx+i) = exp(-0.5 * N_sq *X(i)) * (N_sq / (2 * sigma)* sin(m*X(i)) + m / sigma * cos(m *X(i)) )*cos(sigma*0.125); % Density
    
    % Initial Energy
    H0 = H0 + 0.5*(dx/R_0_centre(i))*(U(i).*U(i) + U(Nx+i).*U(Nx+i));
end 


%% Build flux matrices
% First we construct the system d/dt(U) = DIV(U),
% DIV(U) will be the discrete divergence .
% In this example DIV is the sum of the numerical line fluxes for
% each element.

% Internal cell contributions - this will form a tridiagonal system
for i = 2:Nx-1
    % Flux for velocity U, which depends on R
    DIV(i, Nx + i -1) =  theta / ( dx * R_0_centre(i -1 )) * R_0_face( i ) ;
    DIV(i, Nx + i ) =  -theta / ( dx * R_0_centre(i ))* R_0_face( i ) ...
                   + (1- theta)/ ( dx * R_0_centre(i  ))* R_0_face( i+1 ) ...
                   - dR_0_centre(i)/R_0_centre(i);
    DIV(i, Nx + i +1) =  -(1 - theta)/ ( dx * R_0_centre(i + 1 )) * R_0_face( i + 1 ) ; 
    % Flux for density R, which depends on U
    DIV(Nx + i,  i -1) = (1 - theta)/ ( dx * R_0_centre(i -1 )) * R_0_face( i ) ;
    DIV(Nx + i,  i ) =  theta/ ( dx * R_0_centre(i  ))* R_0_face( i ) ...
                  - (1-theta)/ ( dx * R_0_centre(i ))* R_0_face( i+1 ) ...
                  +  dR_0_centre(i)/R_0_centre(i);
    DIV(Nx + i,  i +1) = - theta/ ( dx * R_0_centre(i + 1 )) * R_0_face( i + 1 ) ; 
end    

% Solid wall boundaries
% Left wall
DIV(1,Nx+1) =  (1-theta)/ ( dx * R_0_centre(1 )) * R_0_face( 2 ) -  dR_0_centre(1)/R_0_centre(1);
DIV(1,Nx+2) = - (1-theta)/ ( dx * R_0_centre(2 )) * R_0_face( 2 );
DIV(Nx+1,1) = -(1-theta)/ ( dx * R_0_centre(1 )) * R_0_face( 2 ) +  dR_0_centre(1)/R_0_centre(1);
DIV(Nx+1,2) = -theta/ ( dx * R_0_centre(2 )) * R_0_face( 2 );
% Right wall
DIV(Nx,2*Nx-1) = theta/ ( dx * R_0_centre(Nx -1  )) * R_0_face( Nx );
DIV(Nx,2*Nx) = -theta/ ( dx * R_0_centre(Nx )) * R_0_face( Nx ) -  dR_0_centre(Nx)/R_0_centre(Nx);
DIV(2*Nx,Nx-1) = (1-theta)/ ( dx * R_0_centre(Nx -1  )) * R_0_face( Nx );
DIV(2*Nx,Nx) = theta/ ( dx * R_0_centre(Nx )) * R_0_face( Nx ) +  dR_0_centre(Nx)/R_0_centre(Nx);

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
Inverse = P\Q;

%F(ceil(periods/dt)) = struct('cdata',[],'colormap',[]);
%% Timestep Loop
count_energy = 0;
figure
while t < periods
    t =t+dt;
    H_local = 0;
    % Calculate  U^(n+1) from  U^n and advance in time
    U = Inverse*U;
    % Calculate Energy
    count_energy = count_energy+1;
    for i = 1:Nx
        H_local = H_local + 0.5*(dx/R_0_centre(i))*(U(i).*U(i) + U(Nx+i).*U(Nx+i));
    end
    H(count_energy,1) = H_local;
    
    % Plot velocity 
    plot(X,U(1:Nx))
    title(['Solution after ', num2str(dt*count_energy), ' periods - Velocity u'])
    axis([  a L -1 1])
    xlabel('x')
    ylabel('u')
    pause(0.01)
    hold off
    
    
    
end
%F(count_energy)=getframe(gcf);
%movie2avi(F, 'velocity.avi')

% Plot Error in Energy
figure
hold all
set(gca,'fontsize', 12)
plot(linspace(0,periods,periods/dt), (H-H0)/H0)
title('Relative Error in Energy')
xlabel('Number of Periods')
hold off

