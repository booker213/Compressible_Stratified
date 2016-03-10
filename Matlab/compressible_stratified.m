clc  
clear all 
close all


%% mesh constants
a = 0; % Start of domain
Lx = 1; % Length of domain in x direction
Lz = 1; % Length of domain in z direction
Nx=4; % Number of elements in x direction
Nz=4; % Number of elements in z direction
k= Nx*Nz; % Number of elements in domain
k4 = 4*k; % Number of local discrete variables in domain
N = 2; % Buoyancy frequency
g = 1; % gravity constant
c0 =  1; % Speed of sound constant 
sigma = sqrt( 0.5* ( 8*pi^2 + 0.25* (N+1)^2 + sqrt((8*pi^2 + 0.25* (N+1)^2)^2-16*pi^2*N^2))); % initial condition - wave frequency
period = 1/sigma; % initial condition - wave period
dt= 1/16; % time step
%dt = period / (10*Nx^2); % time step
t=0; % initial time
tend = 10*period; % end time
theta = 0.5; % flux constant, 0 < theta < 1
steps = tend/dt; % number of time steps taken
% Storage for  discrete energy
energy= zeros(ceil(steps+1),0);

%% Storage for variables 
% Storage for solution vector
% U = [U,W,R,P]^T
U = zeros(k4,0);
% Storage for flux matrices
% S contains the source terms for each element in the domain
% F contains the flux contributions for each element in the domain
% DIV - the discrete divergence matrix is created by
% DIV = S + F
S = zeros(k4, k4 );
F = zeros(k4, k4 );
% Storage for initial condition.

uu = zeros(Nz, Nx);
uw = zeros(Nz, Nx);
r = zeros(Nz, Nx);
p = zeros(Nz, Nx);
% Storage for local discrete energy
H = zeros(Nz, Nx);
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



%% Background density 
r_0 = exp(-3.*Z);

% The flux for the pressure and density
% contain a vertical velocity term that has a dependance on the 
% ratio of the background density in either cell.
% As our chosen stratification is uniform we will calculate these ratios
% here.
r_0down  = r_0(1,1)/r_0(2,1) ;
r_0up  = r_0(2,1)/r_0(1,1) ;


% matlab matrix (columns, rows)
%% Initial condition
for i= 1: Nz
   for j = 1:Nx
uu(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi*Z(i,j)) - (N-1)*0.5*sin(2*pi*Z(i,j)))*sin(2*pi*X(i,j))*sin(0.1);
uw(i,j) = exp(-0.5*(N + 1)*Z(i,j))*sin(2*pi*Z(i,j))*cos(2*pi*X(i,j))*sin(0.1);
r(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(0.1);
p(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(0.1);

   end
end

% Evaluate initial discrete energy    
for i= 1: Nz
  for j = 1:Nx
H(i,j) = dv*((0.5/r_0(i,j))*(uu(i,j)^2+uw(i,j)^2) + (0.5/(r_0(i,j)* N))*(r(i,j) - p(i,j))^2 + 0.5* p(i,j)^2/r_0(i,j));
  end
end
% Assign initial discrete energy to storage
energy(1,1) = sum(sum(H));

%% Initial Condition
% Reshape initial condition into column vector
% U = [U,W,R,P]^T
u_ = reshape(uu', [k, 1]);
w_ = reshape(uw', [k, 1]);
r_ = reshape(r', [k, 1]);
p_ = reshape(p', [k, 1]);
U = vertcat(u_,w_,r_,p_);

%% Make Source matrix S
% W source loop
% W has source terms depending on R and P
count_source_w = 2*k;
% Density source for W
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = (g^3/(N*c0) - 3*(g^2/(N)));
end
% Pressure source for W
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = (-g^3/(N*c0^2) - g/c0 + 3*g^2/(N/c0) + 3);

end

% P source loop
% P has source terms depending on W
count_source_p = k;
for j = 3*k+1:4*k
count_source_p = count_source_p +1;
S(j,count_source_p) = (g - 3*c0);

end

%% Make flux matrix F
% reminder : 
% for an entry in the matrix F(i,j)
% Relating that term to its neighbours are the following
%  +/- dx -> +/- 1
% +/-dz -> +,- Nx


% Flux for U
% Flux only contains the horizontal face contributions of the pressure P
count_col=3*k; % acts on pressure
count_row=0; 
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left corner
                F(count_row, count_col) = (1-theta) ;
                F(count_row, count_col+1) = -(1-theta) ;                
           % bottom right corner
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col-1) = (theta) ;
               F(count_row, count_col) = -theta ;   
           % bottom row - z = 0
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = theta  ;
               F(count_row, count_col) = (1-theta) - theta;
               F(count_row, count_col+1) =  -(1-theta);   
           elseif (i== Nz ) && (j==1)
           % top left corner
                F(count_row, count_col) = (1-theta) ;
                F(count_row, count_col+1) = -(1-theta) ;   
           % top right corner
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-1) = (theta) ;
               F(count_row, count_col) = -theta  ;
              
           % top row - z = Lz
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = theta ;
               F(count_row, count_col) = (1-theta) - theta ;
               F(count_row, count_col+1) = -(1-theta) ;  
            % left wall - x = 0
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col) = (1-theta) ;
                F(count_row, count_col+1) = -(1-theta) ;   
             % right wall - x = Lx
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               
               F(count_row, count_col-1) = (theta) ;
               F(count_row, count_col) = -theta  ;
           else % internal elements
               F(count_row, count_col-1) = theta ;
               F(count_row, count_col) = (1-theta) - theta  ;
               F(count_row, count_col+1) = -(1-theta);   
           
           
       end
   end
end


% Flux for W
% Flux only contains the vertical face contributions of the pressure P
count_col=3*k; % acts on pressure
count_row=k;
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col) = (1-theta);
                F(count_row, count_col+Nx) = -(1-theta) ;                
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
                F(count_row, count_col) = (1-theta) ;
                F(count_row, count_col+Nx) = -(1-theta) ;  
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
                F(count_row, count_col) = (1-theta) ;
                F(count_row, count_col+Nx) = -(1-theta) ;   
           elseif (i== Nz ) && (j==1)
           % top left
               F(count_row, count_col-Nx) = (theta) ;
               F(count_row, count_col) = -theta  ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-Nx) = (theta) ;
               F(count_row, count_col) = -theta ;
              
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-Nx) = (theta) ;
               F(count_row, count_col) = -theta  ;
               
           % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = theta  ;
               F(count_row, count_col) = (1-theta) - theta ;
               F(count_row, count_col+Nx) = -(1-theta) ;     
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = theta ;
               F(count_row, count_col) = (1-theta) - theta;
               F(count_row, count_col+Nx) = - (1-theta) ;  
           else
               F(count_row, count_col-Nx) = theta ;
               F(count_row, count_col) = (1-theta) - theta  ;
               F(count_row, count_col+Nx) = -(1-theta) ;   
           
           
       end
   end
end

% Flux for R
% Flux  contains the vertical face contributions of  W
% and horizontal face contributions of U
count_col_u=0;
count_col_w=k;
count_row=2*k;
for i = 1:Nz
   for j= 1:Nx
        count_col_u = count_col_u+1; % terms for U 
        count_col_w = count_col_w+1; % terms for W
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta ;  
                
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta  ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta;     
               
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta ;  
               
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta ; 
                
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta;%*r_0down ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta ; 
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = (theta);%*r_0down ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta  ;  
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta;%*r_0down ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta  ;
                
                F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta  ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta ;  
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta  ;   
           else
               F(count_row, count_col_u-1) = 1 - theta  ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta  ;   
           
               F(count_row, count_col_w-Nx) = 1 - theta  ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta ;   
           
           
       end
   end
end

% Flux for P
% Flux  contains the vertical face contributions of  W
% and horizontal face contributions of U
count_col_u=0;
count_col_w=k;
count_row=3*k;
for i = 1:Nz
   for j= 1:Nx
        count_col_u = count_col_u+1; % terms for U 
        count_col_w = count_col_w+1; % terms for W
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta ;  
                
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta  ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta;     
               
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta ;  
               
                F(count_row, count_col_w) = -(1-theta);%*r_0up ;
                F(count_row, count_col_w+Nx) = -theta ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta ; 
                
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta;%*r_0down ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta ; 
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = (theta);%*r_0down ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta  ;  
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta;%*r_0down ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = -(1-theta) ;
                F(count_row, count_col_u+1) = -theta  ;
                
                F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta  ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = (1-theta) ;
               F(count_row, count_col_u) = theta ;  
               
               F(count_row, count_col_w-Nx) = (1-theta) ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta  ;   
           else
               F(count_row, count_col_u-1) = 1 - theta  ;
               F(count_row, count_col_u) = theta - (1-theta) ;
               F(count_row, count_col_u+1) = -theta  ;   
           
               F(count_row, count_col_w-Nx) = 1 - theta  ;
               F(count_row, count_col_w) = theta - (1-theta);%theta*r_0down - (1-theta)*r_0up ;
               F(count_row, count_col_w+Nx) = -theta ;   
           
           
       end
   end
end

%% Make DIV Matrix - discrete divergence
% System so far is of the form
% d/dt U = DIV(U)


DIV = S + F;

%% Make P, Q matrix
% Apply implicit midpoint timestepper
% System to solve is P U^n+1 = Q U^n

P = eye(k4) - dt*0.5*DIV;
Q = eye(k4) + dt*0.5*DIV;
% P, Q are not time dependant so we create their inverse here.
Inverse = (P\Q);

%% Solve loop
count_energy=1;
while t < tend
count_energy = count_energy+1;
% Apply timestepper
U = Inverse*U;
% Update time
t = t +dt;

%% Decompose numerical solution for energy and error analysis
% Return column vector to 4 separate matrices for each variable
% Used to plot and work out discrete energy

count_col =0;

for i= 1: Nz
  for j = 1:Nx
        count_col=count_col+1;
uu(i,j) = U(count_col); 
    end
end
    for i= 1: Nz
        for j = 1:Nx
        count_col=count_col+1;
uw(i,j) = U(count_col); 
    end
end
    for i= 1: Nz
        for j = 1:Nx
        count_col=count_col+1;
r(i,j) = U(count_col); 
    end
end
    for i= 1: Nz
        for j = 1:Nx
        count_col=count_col+1;
p(i,j) = U(count_col); 
    end
    end
% Energy Evaluation
for i= 1: Nz
  for j = 1:Nx
H(i,j) = dv*((0.5/r_0(i,j))*(uu(i,j)^2+uw(i,j)^2) + (0.5/(r_0(i,j)* N))*(r(i,j) - p(i,j))^2 + 0.5* p(i,j)^2/r_0(i,j));
  end
end
energy(count_energy,1) = sum(sum(H));

end

H_0 =  energy(1,1);
 figure 
plot(linspace(0,tend, count_energy),(energy-H_0)/H_0)
title('Relative error in energy')
xlabel('Time')
ylabel('energy')
