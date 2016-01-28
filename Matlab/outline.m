clc 
clear all

% mesh constants
a = 0;
Lx = 1;
Lz = 1;
Nx=5;
Nz=5;
k= Nx*Nz;
k4 = 4*k;
N = 2;
dt= 1/16;
t=0;

theta = 0.5;
%% Random theta
% theta_ = rand(k, 1);
% theta = repmat( theta_, 4);

% Storage for inital density and density gradient
r_000=zeros(Nx, Nz);
r_00 = zeros(k,0);
r_0 = zeros(k4,0);
dr_0=zeros(k4, 0);
% Storage for variables 
U= zeros(k4,0);
S = zeros(k4, k4 );
F = zeros(k4, k4 );
PB = zeros(k4, k4 );
P = zeros(k4, k4 );
Q = zeros(k4, k4 );
uu_initial = zeros(Nx, Nz);
uw_initial = zeros(Nx, Nz);
r_initial = zeros(Nx, Nz);
p_initial = zeros(Nx, Nz);
% mesh information
x_nodes = linspace(a,Lx,Nx+1);
z_nodes = linspace(a,Lz,Nz+1);

dx = x_nodes(2) - x_nodes(1);
dz = z_nodes(2) - z_nodes(1);
dv= dx*dz;

x_centres = a+dx/2 : dx : Lx;
z_centres = a+dz/2 : dz : Lz;

[X,Z] = meshgrid(x_centres, z_centres);

sigma = sqrt( 0.5* ( 8*pi^2 + 0.25* (N+1)^2 + sqrt((8*pi^2 + 0.25* (N+1)^2)^2-16*pi^2*N^2)));

% Background density 
r_000 = exp(-3.*Z);
% 
% count_r = 0;
% for j = 1:Nz
% for i = 1:Nx
% count_r = count_r + 1;
% r_00(count_r) = r_000(i,j);
% end
% end
% r_0 = repmat(r_00,4);
% dr_0 = -3*r_0;

% Exact Solution 
uu_initial = exp(-0.5*(N + 1).*Z)*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi.*Z) - (N-1)*0.5*sin(2*pi.*Z))*sin(2*pi.*X)*sin(sigma*t+0.1);
uw_initial = exp(-0.5*(N + 1).*Z)*sin(2*pi.*Z)*cos(2*pi.*X)*sin(sigma*t+0.1);
r_initial = exp(-0.5*(N + 1).*Z)*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi.*Z)-2*pi*cos(2*pi.*Z))*cos(2*pi.*X)*cos(sigma*t+0.1);
p_initial = exp(-0.5*(N + 1).*Z)*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi.*Z)-2*pi*cos(2*pi.*Z))*cos(2*pi.*X)*cos(sigma*t+0.1);


H = (0.5./r_000).*(uu_initial.^2+uw_initial.^2) + (0.5./(r_000 .* N)).*(r_initial - p_initial).^2 + 0.5.* p_initial.^2./r_000;


energy_initial = trapz(x_centres, trapz(z_centres,H, 2))

% Initial Condition
% count_ = 0;
% % u loop
% for j = 1:Nz
% for i = 1:Nx
% count_ = count_ + 1;
% U(count_) = uu_initial(i,j);
% end
% end
u_ = reshape(uw_initial', k, 1);
w_ = reshape(uw_initial', k, 1);
r_ = reshape(r_initial', k, 1);
p_ = reshape(p_initial', k, 1);
U = vertcat(u_,w_,r_,p_);
% % w loop 
% for j = 1:Nz
% for i = 1:Nx
% count_ = count_ + 1;
% U(count_) = uw_initial(i,j);
% end
% end
% 
% % r loop
% for j = 1:Nz
% for i = 1:Nx
% count_ = count_ + 1;
% U(count_) = r_initial(i,j);
% end
% end
% 
% % p loop
% for j = 1:Nz
% for i = 1:Nx
% count_ = count_ + 1;
% U(count_) = p_initial(i,j);
% end
% end

% Make Source matrix
% W source loop
for j = k+1:2*k
for i = 2*k+1:3*k
S(j,i) = dv^2*(-2/N);
end 
for i = 3*k+1:4*k
S(j,i) = dv^2*(2/N+2);
end
end

% P source loop
for j = 3*k+1:4*k
for i = k+1:2*k
S(j,i) = dv^2*(-2);
end 
end

% Make flux matrix

count_col=3*k;
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
       % for count_row = 1:k
       count_row =1 ;
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col) = 1 ;
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col) = 2 ;
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
                F(count_row, count_col)= 5 ;
           elseif (i== Nz ) && (j==1)
           % top left
                  F(count_row, count_col) = 3 ;
           % top right
           elseif (i== Nz ) && (j==Nx)
                F(count_row, count_col) = 4 ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
                F(count_row, count_col) = 6 ;
           else
                F(count_row, count_col) = 7 ;
          % end
           
       end
   end
end
% Values
% +1
% (1-theta)/dx
% 0 
% (2*theta-1)/dx
% -1
% -theta / dx
% reminder : i -> +,- 1
% j -> +,- Nx


% Make Poisson Matrix

PB = S - F;

% Make P, Q matrix
% System to solve is P U^n+1 = Q U^n

P = eye(k4) - dt*0.5*PB;
Q = eye(k4) + dt*0.5*PB;

%U = Q*U/P;