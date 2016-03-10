
clc  
clear all 
close all

%% mesh constants
a = 0; % Start of domain
Lx = 1; % Length of domain in x direction
Lz = 1; % Length of domain in z direction
Nx=32; % Number of elements in x direction
Nz=32; % Number of elements in z direction
k= Nx*Nz; % Number of elements in domain
k4 = 4*k; % Number of local discrete variables in domain
N = 2; % Buoyancy frequency
g = 1; % gravity constant
c0 =  1; % Speed of sound constant 
sigma = sqrt( 0.5* ( 8*pi^2 + 0.25* (N+1)^2 + sqrt((8*pi^2 + 0.25* (N+1)^2)^2-16*pi^2*N^2))); % initial condition - wave frequency
period = 1/sigma; % initial condition - wave period
dt= 0.01; % time step
%dt = period / (10*Nx^2); % time step
t=0; % initial time
tend = 100*period; % end time

steps = tend/dt; % number of time steps taken



% Storage for initial condition.

uu = zeros(Nz, Nx);
uw = zeros(Nz, Nx);
r = zeros(Nz, Nx);
p = zeros(Nz, Nx);

% mesh information
x_nodes = linspace(a,Lx,Nx+1);
z_nodes = linspace(a,Lz,Nz+1);
% element widths
dx = x_nodes(2) - x_nodes(1);
dz = z_nodes(2) - z_nodes(1);
% element volume
dv= dx*dz;
% Location of element centres
x_centres = a+dx/2 : dx : Lx;
z_centres = a+dz/2 : dz : Lz;

[X,Z] = meshgrid(x_centres, z_centres);
figure
while t < tend
    

for i= 1: Nz
    for j = 1:Nx
uu(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi*Z(i,j)) - (N-1)*0.5*sin(2*pi*Z(i,j)))*sin(2*pi*X(i,j))*sin(sigma*t+0.1);
uw(i,j) = exp(-0.5*(N + 1)*Z(i,j))*sin(2*pi*Z(i,j))*cos(2*pi*X(i,j))*sin(sigma*t+0.1);
r(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);
p(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);
    end
end

% % Plot exact solution

contourf(X,Z,uw)
title([' Exact Soln of pressure at t = ', num2str(t)])
xlabel(' x ')
ylabel( ' z ')
colorbar
caxis([-1,1])
hold off
pause(0.1)

t = t +dt;

end