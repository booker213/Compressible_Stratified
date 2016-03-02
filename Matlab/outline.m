clc 
clear all
close all
format long
%% mesh constants
a = 0;
Lx = 1;
Lz = 1;
Nx=4;
Nz=4;
k= Nx*Nz;
k4 = 4*k;
N = 2;
sigma = sqrt( 0.5* ( 8*pi^2 + 0.25* (N+1)^2 + sqrt((8*pi^2 + 0.25* (N+1)^2)^2-16*pi^2*N^2)))
period = 1/sigma;
 %dt= 1/16;
dt = period / (10*Nx^2);
t=0;
tend = period;
theta = 0.5;

steps = tend/dt;

energy= zeros(ceil(steps+1),0);
energy_sim = zeros(ceil(steps+1),0);
%% Random theta
% theta_ = rand(k, 1);
% theta = repmat( theta_, 4);

%% Storage for inital density and density gradient
r_000=zeros(Nz, Nx);

%% Storage for variables 
U= zeros(k4,0);
S = zeros(k4, k4 );
F = zeros(k4, k4 );

uu_initial = zeros(Nz, Nx);
uw_initial = zeros(Nz, Nx);
r_initial = zeros(Nz, Nx);
p_initial = zeros(Nz, Nx);
uu = zeros(Nz, Nx);
uw = zeros(Nz, Nx);
r = zeros(Nz, Nx);
p = zeros(Nz, Nx);
error_uu = zeros(Nz, Nx);
error_uw = zeros(Nz, Nx);
error_r = zeros(Nz, Nx);
error_p = zeros(Nz, Nx);
% mesh information
x_nodes = linspace(a,Lx,Nx+1);
z_nodes = linspace(a,Lz,Nz+1);

dx = x_nodes(2) - x_nodes(1);
dz = z_nodes(2) - z_nodes(1);
dv= dx*dz;

x_centres = a+dx/2 : dx : Lx;
z_centres = a+dz/2 : dz : Lz;

[X,Z] = meshgrid(x_centres, z_centres);



%% Background density 
r_000 = exp(-3.*Z);

% matlab matrix (columns, rows)
%% Exact Solution 
    for i= 1: Nz
        for j = 1:Nx
uu_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi*Z(i,j)) - (N-1)*0.5*sin(2*pi*Z(i,j)))*sin(2*pi*X(i,j))*sin(sigma*t+0.1);
uw_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*sin(2*pi*Z(i,j))*cos(2*pi*X(i,j))*sin(sigma*t+0.1);
r_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);
p_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);

    end
end
    for i= 1: Nz
        for j = 1:Nx
H(i,j) = dv*(0.5./r_000(i,j))*(uu_initial(i,j)^2+uw_initial(i,j)^2) + dv*(0.5./(r_000 (i,j)* N)).*(r_initial(i,j) - p_initial(i,j))^2 + dv*0.5.* p_initial(i,j)^2/r_000(i,j);
    end
end
energy(1,1) = sum(sum(H));
H_0 = energy(1,1);
energy_sim(1,1)= energy(1,1);



%% Initial Condition
% Reshape initial condition into column vector
u_ = reshape(uu_initial', [k, 1]);
w_ = reshape(uw_initial', [k, 1]);
r_ = reshape(r_initial', [k, 1]);
p_ = reshape(p_initial', [k, 1]);
U = vertcat(u_,w_,r_,p_);


%% Make Source matrix
% W source loop
count_source_w = 2*k;
% Density source
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = (-2/N);
end
% Pressure source
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = (2/N+2);

end

% P source loop
% Vertical velocity source
count_source_p = k;
for j = 3*k+1:4*k
count_source_p = count_source_p +1;
S(j,count_source_p) = (-2);

end

%% Make flux matrix
% Flux for u
count_col=3*k;
count_row=0;
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col) = ( (theta))/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;                
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = ( - (1-theta))/dx ;   
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta - (1-theta))/dx ;
               F(count_row, count_col+1) = (1-theta)/dx ;   
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col) = (theta )/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;   
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = ( - (1-theta))/dx ;
              
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta - (1-theta))/dx ;
               F(count_row, count_col+1) = (1-theta)/dx ;  
             % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col) = (theta )/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = ( - (1-theta))/dx ;
           else
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta - (1-theta))/dx ;
               F(count_row, count_col+1) = (1-theta)/dx ;   
           
           
       end
   end
end

% Flux for w
count_col=3*k;
count_row=k;
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col) = (theta )/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;                
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
                F(count_row, count_col) = ( - (1-theta))/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;  
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
                F(count_row, count_col) = (theta )/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;   
           elseif (i== Nz ) && (j==1)
           % top left
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) =((1-theta) +theta - (1-theta))/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = ( - (1-theta))/dz ;
              
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = ( - (1-theta))/dz ;
               
                            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (theta - (1-theta))/dz ;
               F(count_row, count_col+Nx) = (1-theta)/dz ;     
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (theta - (1-theta))/dz ;
               F(count_row, count_col+Nx) = (1-theta)/dz ;  
           else
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (2*theta-1)/dz ;
               F(count_row, count_col+Nx) = (1-theta)/dz ;   
           
           
       end
   end
end

% Flux for r
count_col_u=0;
count_col_w=k;
count_row=2*k;
for i = 1:Nz
   for j= 1:Nx
        count_col_u = count_col_u+1;
        count_col_w = count_col_w+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col_u) = (1-theta)/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ;  
                
                F(count_row, count_col_w) = (1-theta)/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (- (1-theta) )/dx ;
               F(count_row, count_col_u) = (-theta)/dx ;     
               
                F(count_row, count_col_w) = (1- theta  )/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- (1-theta)  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;  
               
                F(count_row, count_col_w) = (1-theta  )/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = (1-theta  )/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ; 
                
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = -(theta)/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = (-theta)/dx ; 
               
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = (-theta)/dz ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = -(1 - theta  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = (-theta)/dz ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = (1-theta  )/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ;
                
                F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = -(theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = -(1 - theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
           else
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;   
           
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
           
           
       end
   end
end

% Flux for p
count_col_u=0;
count_col_w=k;
count_row=3*k;
for i = 1:Nz
   for j= 1:Nx
        count_col_u = count_col_u+1;
        count_col_w = count_col_w+1;
        count_row = count_row+1;
      
            if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col_u) = (1-theta)/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ;  
                
                F(count_row, count_col_w) = (1-theta)/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (- (1-theta) )/dx ;
               F(count_row, count_col_u) = (-theta)/dx ;     
               
                F(count_row, count_col_w) = (1- theta  )/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- (1-theta)  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;  
               
                F(count_row, count_col_w) = (1-theta  )/dz ;
                F(count_row, count_col_w+Nx) = (theta)/dz ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = (1-theta  )/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ; 
                
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = -(theta)/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = (-theta)/dx ; 
               
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = (-theta)/dz ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = -(1 - theta  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = (-theta)/dz ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = (1-theta  )/dx ;
                F(count_row, count_col_u+1) = (theta)/dx ;
                
                F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = -(theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = -(1 - theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
           else
               F(count_row, count_col_u-1) = -(1- theta  )/dx ;
               F(count_row, count_col_u) = ((1-theta) - theta )/dx ;
               F(count_row, count_col_u+1) = (theta)/dx ;   
           
               F(count_row, count_col_w-Nx) = -(1- theta  )/dz ;
               F(count_row, count_col_w) = ((1-theta) - theta )/dz ;
               F(count_row, count_col_w+Nx) = (theta)/dz ;   
           
           
       end
   end
end

% reminder : i -> +,- 1
% j -> +,- Nx


%% Make DIV Matrix


DIV = S + F;

%% Make P, Q matrix
% System to solve is P U^n+1 = Q U^n

P = eye(k4) - dt*0.5*DIV;
Q = eye(k4) + dt*0.5*DIV;

Inverse = (P\Q);
%% Solve loop
count_energy=1;
while t < tend
 count_energy = count_energy+1;

U = Inverse*U;
t = t +dt;
%% Exact Solution for timestep
    for i= 1: Nz
        for j = 1:Nx
uu_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi*Z(i,j)) - (N-1)*0.5*sin(2*pi*Z(i,j)))*sin(2*pi*X(i,j))*sin(sigma*t+0.1);
uw_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*sin(2*pi*Z(i,j))*cos(2*pi*X(i,j))*sin(sigma*t+0.1);
r_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);
p_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);

    end
end
for j = 1:Nx
    for i= 1: Nz
H(i,j) = dv*(0.5./r_000(i,j))*(uu_initial(i,j)^2+uw_initial(i,j)^2) + dv*(0.5./(r_000 (i,j)* N)).*(r_initial(i,j) - p_initial(i,j))^2 + dv*0.5.* p_initial(i,j)^2/r_000(i,j);
    end
end
energy(count_energy,1) = sum(sum(H));


%% Decompose numerical solution for energy and error analysis
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
for j = 1:Nx
    for i= 1: Nz
H(i,j) = dv*(0.5/r_000(i,j))*(uu(i,j)^2+uw(i,j)^2) + dv*(0.5/(r_000 (i,j)* N)).*(r(i,j) - p(i,j))^2 + dv*0.5* p(i,j)^2/r_000(i,j);
    end

end
energy_sim(count_energy,1) = sum(sum(H));

end

% % Plot exact solution
% figure
% contourf(X,Z,p_initial)
% title([' Exact Soln of pressure at t = ', num2str(t)])
% xlabel(' x ')
% ylabel( ' z ')
% colorbar
% 
% 
%  figure
% contourf(X,Z,p)
% pause(0.1)
% title([' Numerical Soln of pressure at t = ', num2str(t)])
% colorbar
% xlabel(' x ')
% ylabel( ' z ')


 error_uu= uu - uu_initial;
 error_uw= uw - uw_initial;
 error_r= r - r_initial;
 error_p= p - p_initial;
 
 % L_infty error
 uu_infty = max(abs(error_uu(:)));
 uw_infty = max(abs(error_uw(:)));
 r_infty = max(abs(error_r(:)));
 p_infty = max(abs(error_p(:)));
 
  % L_2error
 uu_2 = sqrt (sum(sum (error_uu.*error_uu)));
 uw_2 = sqrt (sum(sum (error_uw.*error_uw)));
 r_2 = sqrt (sum (sum(error_r.*error_r)));
 p_2 = sqrt (sum( sum(error_p.*error_p)));

 
 figure 
plot(linspace(0,tend, count_energy),(energy_sim-H_0)/H_0)
title('Relative error in energy')
xlabel('Periods')
ylabel('energy')
