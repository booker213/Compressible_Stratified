clc 
clear all
format long
% mesh constants
a = 0;
Lx = 1;
Lz = 1;
Nx=3;
Nz=4;
k= Nx*Nz;
k4 = 4*k;
N = 2;
dt= 1/16;
t=0;
tend = t+2*dt;
theta = 0.4;
%% Random theta
% theta_ = rand(k, 1);
% theta = repmat( theta_, 4);

% Storage for inital density and density gradient
r_000=zeros(Nz, Nx);
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
uu_initial = zeros(Nz, Nx);
uw_initial = zeros(Nz, Nx);
r_initial = zeros(Nz, Nx);
p_initial = zeros(Nz, Nx);
% mesh information
x_nodes = linspace(a,Lx,Nx+1);
z_nodes = linspace(a,Lz,Nz+1);

dx = x_nodes(2) - x_nodes(1);
dz = z_nodes(2) - z_nodes(1);
dv= dx*dz;

x_centres = a+dx/2 : dx : Lx;
z_centres = a+dz/2 : dz : Lz;

[X,Z] = meshgrid(x_centres, z_centres);

sigma = sqrt( 0.5* ( 8*pi^2 + 0.25* (N+1)^2 + sqrt((8*pi^2 + 0.25* (N+1)^2)^2-16*pi^2*N^2)))

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
% matlab matrix (columns, rows)
% Exact Solution 
for j = 1:Nx
    for i= 1: Nz
uu_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(2*pi/(4*pi^2 - sigma^2))*( -2*pi*cos(2*pi*Z(i,j)) - (N-1)*0.5*sin(2*pi*Z(i,j)))*sin(2*pi*X(i,j))*sin(sigma*t+0.1);
uw_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*sin(2*pi*Z(i,j))*cos(2*pi*X(i,j))*sin(sigma*t+0.1);
r_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*((0.5*(N+1)-4*pi^2*N/sigma^2)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);
p_initial(i,j) = exp(-0.5*(N + 1)*Z(i,j))*(sigma/(4*pi^2-sigma^2))*(-0.5*(N-1)*sin(2*pi*Z(i,j))-2*pi*cos(2*pi*Z(i,j)))*cos(2*pi*X(i,j))*cos(sigma*t+0.1);

    end
end
for j = 1:Nx
    for i= 1: Nz
H(i,j) = (0.5./r_000(i,j))*(uu_initial(i,j)^2+uw_initial(i,j)^2) + (0.5./(r_000 (i,j)* N)).*(r_initial(i,j) - p_initial(i,j))^2 + 0.5.* p_initial(i,j)^2/r_000(i,j);
    end
end
energy_initial = trapz(z_centres, trapz(x_centres,H, 2))

% Initial energy only close when n =64^2
% However if we measure drift this should be fine?
% Otherwise more accurate integration scheme will be required.



% Initial Condition
% count_ = 0;
% % u loop
% for j = 1:Nz
% for i = 1:Nx
% count_ = count_ + 1;
% U(count_) = uu_initial(i,j);
% end
% end
u_ = reshape(uu_initial', [k, 1]);
w_ = reshape(uw_initial', [k, 1]);
r_ = reshape(r_initial', [k, 1]);
p_ = reshape(p_initial', [k, 1]);
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
count_source_w = 2*k;
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = dv^2*(-2/N);
end
for j = k+1:2*k
count_source_w = count_source_w+1;
S(j,count_source_w) = dv^2*(2/N+2);

end

% P source loop
count_source_p = k;
for j = 3*k+1:4*k
count_source_p = count_source_p +1;
S(j,count_source_p) = dv^2*(-2);

end

% Make flux matrix
% Flux for u
count_col=3*k;
count_row=0;
for i = 1:Nz
   for j= 1:Nx
        count_col = count_col+1;
        count_row = count_row+1;
      
           if (i== 1 ) && (   j==1)
           % bottom left
                F(count_row, count_col) = (theta -1 )/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;                
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta)/dx ;   
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (2*theta-1)/dx ;
               F(count_row, count_col+1) = (1-theta)/dx ;   
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col) = (theta -1 )/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;   
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta)/dx ;
              
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (2*theta-1)/dx ;
               F(count_row, count_col+1) = (1-theta)/dx ;  
             % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                          F(count_row, count_col) = (theta -1 )/dx ;
                F(count_row, count_col+1) = (1-theta)/dx ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               
                             F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (theta)/dx ;
           else
               F(count_row, count_col-1) = (- theta  )/dx ;
               F(count_row, count_col) = (2*theta-1)/dx ;
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
                F(count_row, count_col) = (theta -1 )/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;                
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
                F(count_row, count_col) = (theta -1 )/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;  
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
                F(count_row, count_col) = (theta -1 )/dz ;
                F(count_row, count_col+Nx) = (1-theta)/dz ;   
           elseif (i== Nz ) && (j==1)
           % top left
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (theta)/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (theta)/dz ;
              
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (theta)/dz ;
               
                            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (2*theta-1)/dz ;
               F(count_row, count_col+Nx) = (1-theta)/dz ;     
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col-Nx) = (- theta  )/dz ;
               F(count_row, count_col) = (2*theta-1)/dz ;
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
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ;  
                
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ;     
               
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;  
               
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ; 
                
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ; 
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ;
                
                F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
           else
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;   
           
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
           
           
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
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ;  
                
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ;  
           % bottom right
           elseif (i== 1 ) &&  (j==Nx)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ;     
               
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ; 
           % bottom row
           elseif  (i== 1 ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;  
               
                F(count_row, count_col_w) = (theta -1 )/dz ;
                F(count_row, count_col_w+Nx) = (1-theta)/dz ;  
           elseif (i== Nz ) && (j==1)
           % top left
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ; 
                
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
           % top right
           elseif (i== Nz ) && (j==Nx)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ; 
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
           % top row
           elseif  (i== Nz ) && (   j~=1) && (j~=Nx)        
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (theta)/dz ;
            % left wall
           elseif (j== 1 ) && (   i~=1) && (i~=Nz)
                F(count_row, count_col_u) = (theta -1 )/dx ;
                F(count_row, count_col_u+1) = (1-theta)/dx ;
                
                F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
             % right wall
           elseif (j== Nx ) && (   i~=1) && (i~=Nz)
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (theta)/dx ;  
               
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
           else
               F(count_row, count_col_u-1) = (- theta  )/dx ;
               F(count_row, count_col_u) = (2*theta-1)/dx ;
               F(count_row, count_col_u+1) = (1-theta)/dx ;   
           
               F(count_row, count_col_w-Nx) = (- theta  )/dz ;
               F(count_row, count_col_w) = (2*theta-1)/dz ;
               F(count_row, count_col_w+Nx) = (1-theta)/dz ;   
           
           
       end
   end
end

% reminder : i -> +,- 1
% j -> +,- Nx


% Make Poisson Matrix

PB = S - F;

% Make P, Q matrix
% System to solve is P U^n+1 = Q U^n

P = eye(k4) - dt*0.5*PB;
Q = eye(k4) + dt*0.5*PB;
 while t < tend
     t = t +dt;
U = inv(P)*Q*U;

 end