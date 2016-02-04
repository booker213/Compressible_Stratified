from firedrake import *

#Create Mesh
m = 16
mesh = UnitSquareMesh(m, m,  quadrilateral=quadrilateral)
order_basis = 0 

# Note x[0] = z
#      x[1] = x
# Attempt to make an easy extension to 3D ?

# Declare timestep 
# Currently replicating Sanders' work
dt = Constant (1/16)
t = 0.0
end_time  = Constant(1/16)

#Define Constants
c_0 = Constant(1.0) # Speed of sound
g = Constant(9.81) # Gravity
N = Constant(2.0) # Bouyancy frequency
theta = Constant(0.5) # Alternating Flux
#Future make theta differ in elements, hopefully with a random seeding
#Theta \in [ 0.1 , 0.9 ] would be a preferable range.


#Define Function Spaces

V = VectorFunctionSpace(mesh, "DG", order_basis) # Dim of VFS = dim mesh, unless specified otherwise.
R = FunctionSpace(mesh, "DG", order_basis)
P = FunctionSpace(mesh, "DG", order_basis)
W =  V*R*P

# Define Background Density
r_0_expression = Expression( "exp(-3.0*x[0])" )
#r_0_expression = Expression( "1.0" )
# Project onto Finite element space
# so we can use it in variational forms
r_0 = Function(R)
r_0.interpolate(r_0_expression)

#Define trial and test functions
(u,r,p) = TrialFunctions(W) # Velocity, density, pressure respectively
(phi, xi, sigma) = TestFunctions(W) # Test functions for velocity, density, pressure respectively

#Function Space for initial conditions
w_ = Function(W)
(u_, r_, p_) = w_.split()

#Set up element normal
n = FacetNormal(mesh)


#Set up boundary conditions
#no normal flow
# ids 1,2 correspond to z=0, z=1
# reminder that cartesians are rotated 90 deg
#           x=1
#     |-------------|
#  z  |             | z
#  =  |             | = 
#  0  |             | 1
#     |             |
#     |-------------|
#          x = 0
bcu_z_1 = DirichletBC(W.sub(0).sub(0), Constant(0.0), 1, method="geometric")
bcu_z_2 = DirichletBC(W.sub(0).sub(0), Constant(0.0), 2, method="geometric")

# ids 3,4 correspond to x=0, x=1

bcu_x_1 = DirichletBC(W.sub(0).sub(1), Constant(0.0), 3, method="geometric")
bcu_x_2 = DirichletBC(W.sub(0).sub(1), Constant(0.0), 4, method="geometric")



#Define Poisson Bracket
#Possibly as a python definition
def Poisson_Bracket(velocity, density, pressure):
   print "blank atm"



# Define exact solution
#omega = Expression(sqrt(0.5*(8*pow(pi,2)+0.25*pow((N+1),2) + sqrt( pow((8*(pi,2) + 0.25*pow((N+1),2)),2) - 16* (pi,2) *N))))
#Dispersion relation fix later
#Current approach does not work
#Externally determine dispersion relation and insert as constant may be appropiate for now
# Matlab value
omega =  8.900185996715988

#Energy is not as expected.
# Currently need rechecking.
#Interpolating now works

#exact_u = Expression( " exp( -0.5*( N + 1 )*x[0] )* ( 2 * pi/( 4 * pi*pi - omega*omega)) * ( -2*pi*cos(2*pi*x[0]) - 0.5*(N-1)*sin(2*pi*x[0]) )* sin(2*pi*x[1])*sin (omega  + 0.1)", N = 2, omega =  8.900185996715988 )

#exact_w =  Expression( " exp( -0.5*( N + 1 )*x[0] )* sin(2*pi*x[0])*cos(2*pi*x[1])*sin(omega  + 0.1)", N = 2, omega =  8.900185996715988  )

exact_r = Expression( " exp( -0.5*( N + 1 )*x[0] )*(omega / ( 4*pi*pi - omega*omega))*( ( 0.5*(N + 1) - (4*pi*pi*N/(omega*omega)))*sin(2*pi*x[0]) - 2*pi*cos(2*pi*x[0]))*cos(2*pi*x[1])*cos(omega + 0.1) " , N = 2, omega =  8.900185996715988)

exact_p = Expression( " exp( -0.5*( N + 1 )*x[0] )*(omega / ( 4*pi*pi - omega*omega))*( -0.5*(N-1)*sin(2*pi*x[0]) - 2*pi *cos(2* pi *x[0]))*cos(2*pi*x[1])*cos(omega + 0.1) ", N = 2, omega =  8.900185996715988)



#Define initial conditions

r_.interpolate(exact_r)
p_.interpolate(exact_p)
#u_.sub(0).interpolate(exact_w)
#u_.sub(1).interpolate(exact_u)
u_.interpolate(Expression([" exp( -0.5*( N + 1 )*x[0] )* sin(2*pi*x[0])*cos(2*pi*x[1])*sin(omega  + 0.1)",  " exp( -0.5*( N + 1 )*x[0] )* ( 2 * pi/( 4 * pi*pi - omega*omega)) * ( -2*pi*cos(2*pi*x[0]) - 0.5*(N-1)*sin(2*pi*x[0]) )* sin(2*pi*x[1])*sin (omega  + 0.1)"], N = 2, omega =  8.900185996715988))

## u_.interpolate(Expression([component_0, component_1]))
# Lawrence fix
# Interpolation to indexed pieces of a VFS is not currently implemented.

#Initial Energy
E_0 = assemble( (0.5*inner((u_),(u_))/r_0 + 0.5*pow(g,2)*pow(( r_ - p_/c_0),2)/(r_0*N) + 0.5* pow(p_,2)/(r_0*c_0))*dx )


# Create linear problem
# a =   ( u*phi + r*xi + p* sigma)*dx - 0.5*dt*Poisson_Bracket( u , r , p )
# L =   ( u_*phi + r_*xi + p_* sigma)*dx + 0.5*dt*Poisson_Bracket( u_ , r_ , p_ )



# visualisation files
u_file = File('./Results/u.pvd')
#u_w_file = File('./Results/u_w.pvd')
density_file = File('./Results/density.pvd')
pressure_file = File('./Results/pressure.pvd')

# output initial conditions
# output does not currently work
# issue must be the initial conditions 
density_file << r_
pressure_file << p_
u_file << u_
#u_u_file << u_.sub(1)

out=Function(W)

#while (t <= end):


#solve(a == L, out, bcs=[bcv_x_1,bcv_x_2,bcv_z_1,bcv_z_2])
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
print "t = ",  t ,". " "Energy, e = ", E_0
