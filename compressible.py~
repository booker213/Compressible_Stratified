from firedrake import *

#Create Mesh
m = 16
mesh = UnitSquareMesh(m, m)
order_basis = 0 

# Note x[0] = z
#      x[1] = x
# Attempt to make an easy extension to 3D ?

# Declare timestep 
# Currently replicating Sanders' work
dt = Constant (1/16)
t = 0.0
end = 1

#Define Constants
c_0 = Constant(1) # Speed of sound
g = Constant(9.81) # Gravity
N = Constant(2) # Bouyancy frequency
theta = Constant(0.5) # Alternating Flux
#Future make theta differ in elements, hopefully with a random seeding
#Theta \in [ 0.1 , 0.9 ] would be a preferable range.

# Define Background Density
r_0 = Expression( "exp(-3.0*x[0]" )

#Define Function Spaces

V = VectorFunctionSpace(mesh, "DG", order_basis) # Dim of VFS = dim mesh, unless specified otherwise.
R = FunctionSpace(mesh, "DG", order_basis)
P = FunctionSpace(mesh, "DG", order_basis)
W =  V*R*P

#Define trial and test functions
(u,r,p) = TrialFunctions(W) # Velocity, density, pressure respectively
(phi, xi, sigma) = TestFunctions(W) # Test functions for velocity, density, pressure respectively

#Function Space for initial conditions
w_ = Function(W)
(u_, r_, p_) = w_.split()

#Set up element normal
n = FacetNormal(mesh)


#Set up boundary conditions
for i in range(1,4):
   bcu_z_i = DirichletBC(W.sub(0).sub(0), Constant(0), i, method="geometric")
   bcu_x_i = DirichletBC(W.sub(0).sub(1), Constant(0), i, method="geometric")



#Define Poisson Bracket
#Possibly as a python definition
def Poisson_Bracket(velocity, density, pressure):
   print end



# Define exact solution
#omega = Expression(sqrt(0.5*(8*pow(pi,2)+0.25*pow((N+1),2) + sqrt( pow((8*(pi,2) + 0.25*pow((N+1),2)),2) - 16* (pi,2) *N))))
#Dispersion relation fix later
#Current approach does not work
#Externally determine dispersion relation and insert as constant may be appropiate for now
omega = 1.0
exact_u = Expression( " exp( -0.5*( N + 1 )*x[1] )* ( 2 * pi/( 4 * pi^2 - omega^2)) * ( -2*pi*cos(2*pi*x[0]) - 0.5*(N-1)*sin(2*pi*x[0]) )* sin(2*pi*x[1]*sin (omega * t + 0.1)" , t=t)
exact_w =  Expression( " exp( -0.5*( N + 1 )*x[1] )* sin(2*pi*x[0])*cos(2*pi*x[1])*sin(omega *t + 0.1)" , t=t)
exact_r = Expression( " exp( -0.5*( N + 1 )*x[1] )*(omega / ( 4*pi^2 - omega^2))*( ( 0.5*(N + 1) - (4*pi^2*N/omega^2))*sin(2*pi*x[0]) - 2*pi*cos(2*pi*x[0]))*cos(2*pi*x[1])*cos(omega * t +0.1) " , t=t)
exact_p = Expression( " exp( -0.5*( N + 1 )*x[1] )*(omega / ( 4*pi^2 - omega^2))*( -0.5*(N-1)*sin(2*pi*x[0]) - 2*pi *cos(2* pi *x[0]))cos(2*pi*x[1])*cos(omega * t + 0.1) ", t=t)

#Define initial conditions
r_.interpolate(exact_r)
p_.interpolate(exact_p)
u_.sub(0).interpolate(exact_w)
u_.sub(1).interpolate(exact_u)



# Create linear problem
# a =   ( u*phi + r*xi + p* sigma)*dx - 0.5*dt*Poisson_Bracket( u , r , p )
# L =   ( u_*phi + r_*xi + p_* sigma)*dx + 0.5*dt*Poisson_Bracket( u_ , r_ , p_ )



# visualisation files
u_u_file = File('./Results/u_u.pvd')
u_w_file = File('./Results/u_w.pvd')
density_file = File('./Results/density.pvd')
pressure_file = File('./Results/pressure.pvd')

out=Function(W)

#while (t <= end):


#solve(a == L, out, bcs=[bcv_x_1,bcv_x_2,bcv_x_3,bcv_x_4,bcv_z_1,bcv_z_2,bcv_z_3,bcv_z_4])
#velocity, density, pressure = out.split( )

#density_file << density
#pressure_file << pressure
#u_w_file << velocity.sub(0)
#u_u_file << velocity.sub(1)

#Assemble Energy
#E = assemble( (0.5*(velocity)**2/r_0 + 0.5*g*( density - pressure/c_0)**2/(r_0*N) + 0.5* pressure**2/(r_0*c_0))*dx )
#t+= dt



#Compile test
#remove later on
print t,"Seems to have worked"
