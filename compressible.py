from firedrake import *

#Create Mesh
m=16
mesh = UnitSquareMesh(m, m)
order_basis = 0 

# Note x[0] = x
#      x[1] = z

#Declare timestep 
timestep = 1/16
t=0.0
end=1
#Define Constants
c_0 = 1 # speed of sound
g = 9.81 # gravity
N = 2 # Bouyancy frequency
pi= 3.41
# Define Background Density
r_0 = Expression( "exp(-3.0*x[1]" )

#Define Function Spaces

V = VectorFunctionSpace(mesh, "DG", order_basis) # Dim of VFS = dim mesh, unless specified otherwise.
R = FunctionSpace(mesh, "DG", order_basis)
P = FunctionSpace(mesh, "DG", order_basis)
W =  V*R*P

#Define trial and test functions
(u,r,p) = TrialFunctions(W) # Velocity, density, pressure respectively
(phi, xi, sigma) = TestFunctions(W)


#Set up element normal
n = FacetNormal(mesh)
#un = 0.5*(dot(u0, n) + abs(dot(u0, n)))

#Set up boundary conditions
for i in range(1,4):
   bcv_x_i = DirichletBC(W.sub(0).sub(0), Constant(0), i, method="geometric")
   bcv_z_i = DirichletBC(W.sub(0).sub(1), Constant(0), i, method="geometric")

#outfile = File("./Results/compressible.pvd")

#Define Poisson Bracket
#Possibly as a python definition

# Define exact solution
# should really define 2pi as a constant
#omega = Expression(sqrt(0.5*(8*pow(pi,2)+0.25*pow((N+1),2) + sqrt( pow((8*(pi,2) + 0.25*pow((N+1),2)),2) - 16* (pi,2) *N))))
#Dispersion relation fix later
#Current approach does not work
#Externally determine dispersion relation and insert as constant may be appropiate for now
omega = 1.0
exact_u = Expression( " exp( -0.5*( N + 1 )*x[1] )* ( 2 * pi/( 4 * pi^2 - omega^2)) * ( -2*pi*cos(2*pi*x[1]) - 0.5*(N-1)*sin(2*pi*x[1]) )* sin(2*pi*x[0]*sin (omega * t + 0.1)" )
exact_w =  Expression( " exp( -0.5*( N + 1 )*x[1] )* sin(2*pi*x[1])*cos(2*pi*x[0])*sin(omega *t + 0.1)" )
exact_r = Expression( " exp( -0.5*( N + 1 )*x[1] )*(omega / ( 4*pi^2 - omega^2))*( ( 0.5*(N + 1) - (4*pi^2*N/omega^2))*sin(2*pi*x[1]) - 2*pi*cos(2*pi*z))*cos(2*pi*x[0])*cos(omega * t +0.1) " )
exact_p = Expression( " exp( -0.5*( N + 1 )*x[1] )*(omega / ( 4*pi^2 - omega^2))*( -0.5*(N-1)*sin(2*pi*x[1]) - 2*pi *cos(2* pi *x[1]))cos(2*pi*x)*cos(omega * t + 0.1) ")

print t,"Seems to have worked"
# visualisation files
u_file = File('./Results/u.pvd')
w_file = File('./Results/w.pvd')
density_file = File('./Results/density.pvd')
pressure_file = File('./Results/pressure.pvd')
#while (t <= end):

#w=Function(W)
#solve(a == L, w, bcs=[bcv_x_1,bcv_x_2,bcv_x_3,bcv_x_4,bcv_z_1,bcv_z_2,bcv_z_3,bcv_z_4])
#velocity, density, pressure = w.split( )

#Assemble Energy
#E = assemble( (0.5*(velocity)**2/r_0 + 0.5*g*( density - pressure/c_0)**2/(r_0*N) + 0.5* pressure**2/(r_0*c_0))*dx )

