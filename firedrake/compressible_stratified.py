from firedrake import *




#Create Mesh
# Chosen mesh is a squrare with m square elements.
Nx = 16
Nz = 16
mesh = UnitSquareMesh(Nx, Nz,  quadrilateral=quadrilateral)
order_basis = 0 

# Note x[0] = z
#      x[1] = x
# In effect we have stratified the fluid in the x direction 
# instead of the z direction. This is fine as long as gravity 
# acts in the z direction.


# Setup time conditions
t = 0.0
end = 1.

# Declare timestep 
dt = 1./16.

#Define Constants
c_0 = Constant(1.0) # Speed of sound
g = Constant(1.) # Gravity
N = Constant(2.0) # Bouyancy frequency
theta = Constant(0.5) # Alternating Flux
#Future make theta differ in elements, hopefully with a random seeding
#Theta \in [ 0.1 , 0.9 ] would be a preferable range.


#Define Function Spaces
# and create mixed function space
V = VectorFunctionSpace(mesh, "DG", order_basis) # Dim of VFS = dim mesh, unless specified otherwise.
R = FunctionSpace(mesh, "DG", order_basis)
P = FunctionSpace(mesh, "DG", order_basis)
W =  V*R*P

# Define Background Density
r_0 = Function(R)
dr_0 = Function(R)
r_0.interpolate(Expression( "exp(-3.0*x[0])" ))
dr_0.interpolate(Expression( "-3.0*exp(-3.0*x[0])" ))

# Output background density to check orientation

backgrounddensity_file = File('./Results/background_density.pvd')
backgrounddensity_file << r_0

#Define trial and test functions
(u,rho,p) = TrialFunctions(W) # Velocity, density, pressure respectively
(dFdu_vec, dFdrho, dFdp) = TestFunctions(W) # Test functions for velocity, density, pressure respectively

#Function Space for initial conditions

w_ = Function(W)
(u0, rho0, p0) = split(w_)
w = Function(W)


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




# Define exact solution
#omega = Expression(sqrt(0.5*(8*pow(pi,2)+0.25*pow((N+1),2) + sqrt( pow((8*(pi,2) + 0.25*pow((N+1),2)),2) - 16* (pi,2) *N))))
#Dispersion relation fix later
#Current approach does not work
#Externally determine dispersion relation and insert as constant may be appropiate for now

# Matlab value for wave frequency
omega =  8.900185996715988


#Define initial conditions
(u0, rho0, p0) = w_.split()

exact_rho = Expression( " exp( -0.5*( N + 1 )*x[0] )*(omega / ( 4*pi*pi - omega*omega))*( ( 0.5*(N + 1) - (4*pi*pi*N/(omega*omega)))*sin(2*pi*x[0]) - 2*pi*cos(2*pi*x[0]))*cos(2*pi*x[1])*cos(0.1) " , N = 2, omega =  8.900185996715988)

exact_p = Expression( " exp( -0.5*( N + 1 )*x[0] )*(omega / ( 4*pi*pi - omega*omega))*( -0.5*(N-1)*sin(2*pi*x[0]) - 2*pi *cos(2* pi *x[0]))*cos(2*pi*x[1])*cos( 0.1) ", N = 2, omega =  8.900185996715988)

rho0.interpolate(exact_rho)
p0.interpolate(exact_p)
u0.interpolate(Expression([" exp( -0.5*( N + 1 )*x[0] )* sin(2*pi*x[0])*cos(2*pi*x[1])*sin( 0.1)",  " exp( -0.5*( N + 1 )*x[0] )* ( 2 * pi/( 4 * pi*pi - omega*omega)) * ( -2*pi*cos(2*pi*x[0]) - 0.5*(N-1)*sin(2*pi*x[0]) )* sin(2*pi*x[1])*sin ( 0.1)"], N = 2, omega =  8.900185996715988))

## u_.interpolate(Expression([component_0, component_1]))


#Initial Energy
E0 = assemble( ((inner(u0,u0)/r_0 + (rho0**2 - 2*rho0*p0/c_0 + (p0/c_0)**2)/(r_0*N) + (p0**2)/(r_0*c_0) )*0.5)*dx)


# Create linear problem
# a =   ( u*dFdu_vec + r*dFdr + p* dFdp)*dx - 0.5*dt*Poisson_Bracket( u , r , p )
# L =   ( u_*dFdu_vec + r_*dFdr + p_* dFdp)*dx + 0.5*dt*Poisson_Bracket( u_ , r_ , p_ )


#Define varitional derivatives
dHdu = u/r_0
dHdrho = (rho - p/c_0)/(N*r_0)
dHdp = p/(c_0*r_0) + ( p/(c_0**2) - rho/c_0)/(r_0*N)

dHdu0 = u0/r_0
dHdrho0 = (rho0 - p0/c_0)/(N*r_0)
dHdp0 = p0/(c_0*r_0) + ( p0/(c_0**2) - rho0/c_0)/(r_0*N)

#Define Poisson Bracket
#Possibly as a python definition

a0 = (dot(u, dFdu_vec) + rho*dFdrho + p*dFdp)*dx
a1 = (- dot (grad(dHdrho*r_0), dFdu_vec) + dot(grad(r_0*dFdrho), dHdu))*dx
a2 = ( dr_0*(dHdrho*dFdu_vec[0] - dFdrho*dHdu[0]))*dx
a3 = (g*r_0*(dFdp*dHdu[0] - dHdp*dFdu_vec[0] ) )*dx
a4 = (- dot (grad(c_0*r_0*dHdp), dFdu_vec) + dot(grad(c_0*r_0*dFdp), dHdu))*dx
a5 = ( -jump(r_0*dHdrho)*dot((1-theta)*dFdu_vec('-')+ theta*dFdu_vec('+'), n('-')))*dS
a6 = ( jump(r_0*dFdrho)*dot((1-theta)*dHdu('-')+ theta*dHdu('+'), n('-')))*dS
a7 = ( -jump( c_0*r_0*dHdp )*dot((1-theta)*dFdu_vec('-')+ theta*dFdu_vec('+'), n('-')))*dS
a8 = ( jump(c_0*r_0*dFdp)*dot((1-theta)*dHdu('-')+ theta*dHdu('+'), n('-')))*dS

a = a0 - 0.5*dt*( a1 + a2 +  a3 + a4 + a5 + a6 + a7 + a8)

L0 = (dot(u0, dFdu_vec) + rho0*dFdrho + p0*dFdp)*dx
L1 = (- dot (grad(dHdrho0*r_0), dFdu_vec) + dot(grad(r_0*dFdrho), dHdu0))*dx
L2 = ( dr_0*(dHdrho0*dFdu_vec[0] - dFdrho*dHdu0[0]))*dx
L3 = (g*r_0*(dFdp*dHdu0[0] - dHdp0*dFdu_vec[0] ) )*dx
L4 = (- dot (grad(c_0*r_0*dHdp0), dFdu_vec) + dot(grad(c_0*r_0*dFdp), dHdu0))*dx
L5 = ( -jump(r_0*dHdrho0)*dot((1-theta)*dFdu_vec('-')+ theta*dFdu_vec('+'), n('-')))*dS
L6 = ( jump(r_0*dFdrho)*dot((1-theta)*dHdu0('-')+ theta*dHdu0('+'), n('-')))*dS
L7 = ( -jump( c_0*r_0*dHdp0 )*dot((1-theta)*dFdu_vec('-')+ theta*dFdu_vec('+'), n('-')))*dS
L8 = ( jump(c_0*r_0*dFdp)*dot((1-theta)*dHdu0('-')+ theta*dHdu0('+'), n('-')))*dS

L = L0 + 0.5*dt*( L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8)

# Note that we have no boundary surface terms as n.dh/du = n.df/du = 0 at the boundaries


# Storage for visualisation
outfile = File('./Results/compressible_stratified_results.pvd')

u0.rename("Velocity")
rho0.rename("Density")
p0.rename("Pressure")

# Output initial conditions
outfile.write(u0,rho0,p0, time = t)
# File for energy output
E_file = open('./Results/energy.txt', 'w')

out=Function(W)


while (t < end):
 t+=dt
 
 solve(a == L, out)
 u, rho, p = out.split()

 # Assign appropriate name in results file
 u.rename("Velocity")
 rho.rename("Density")
 p.rename("Pressure")
 
 # Output results
 outfile.write(u, rho,p, time =t)
 
 # Assign output as previous timestep for next time update
 u0.assign(u)
 rho0.assign(rho)
 p0.assign(p)
 #Assemble Energy
 E = assemble( ((inner(u,u)/r_0 + (rho**2 - 2*rho*p/c_0 + (p/c_0)**2)/(r_0*N) + (p**2)/(r_0*c_0) )*0.5)*dx )
 E_file.write('%-10s %-10s\n' % (t,abs((E-E0)/E0)))
 # Print time and energy drift, drift should be around machine precision.
 print t, abs((E-E0)/E0)
# Close energy write
E_file.close()

