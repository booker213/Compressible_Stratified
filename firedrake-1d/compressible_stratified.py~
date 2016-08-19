from firedrake import *

# Script to solve compressible stratified waves in upto three dimensions
# c^2_0 and gravity will be considered to be scaled to be equal to 1.

# Test problem will be the one-dimensional waves considered in the MATLAB scripts.

# Create mesh
# Current mesh is a unit line  with Nx elements.

Nx = 64
mesh = UnitIntervalMesh(Nx)

# Declare timestep
dt = 1./( pow(Nx,2))

# Declare initial and end time
# Period of waves considered in test is 1s
# We will consider 3 periods initially
t = 0.
end = 1.

# Declare order of the basis in the elements
# Test problem will consider order 0 - a finite volume scheme
order_basis = 1
quad_degree = 20
# Declare flux indicator function
theta = Constant(0.5)

N_sq = Constant(3.0)
m = 2 * pi
sigma = pow( 0.25 * pow(N_sq, 2) + pow(m, 2),0.5)
#sigma = 1

# Define Background Density

Ro = FunctionSpace(mesh, "DG", order_basis)

r_0 = Function(Ro)
dr_0 = Function(Ro)
r_0.interpolate(Expression( "exp(-3.0*x[0])" ))
dr_0.interpolate(Expression( "-3.0*exp(-3.0*x[0])" ))


# Declare function spaces on the mesh
# and create mixed function space
# for coupled Poisson Bracket.

V = VectorFunctionSpace(mesh, "DG", order_basis)
R = FunctionSpace(mesh, "DG", order_basis)

W = V*R



# Initial conditions

# Function space
w0 = Function(W)


# Interpolate expressions
u0,rho0 = w0.split()

u0.interpolate(Expression("exp(-0.5* N_sq * x[0])* sin ( m * x[0])*sin(sigma*0.125)", N_sq = N_sq, m = m , sigma = sigma ))

rho0.interpolate(Expression("exp(-0.5* N_sq * x[0])*((N_sq / (2 * sigma))* sin ( m * x[0])+ m / sigma * cos ( m * x[0]))*cos(sigma*0.125)", N_sq = N_sq, m = m , sigma = sigma ))

# Assemble initial energy
E0 = assemble ( (0.5/r_0*(inner(u0,u0) + rho0**2))*dx(degree = quad_degree))

print E0

# Set up internal boundary normal on mesh 
# Needed for numerical fluxes in bilinear form

n = FacetNormal(mesh)


# Create bilinear form for linear solver
# Bilinear problem is of the form 
# a(u,v) = L(v)
# Using our Poisson bracket and an implicit midpoint rule
# we see that
# a(u,v) = u^{n+1}*dFdu_vec + rho^{n+1)*dFdrho - 0.5*dt*PB(u^{n+1}, rho^{n+1})
# L(v) = u^{n}*dFdu_vec + rho^{n)*dFdrho + 0.5*dt*PB(u^{n}, rho^{n})


# We note that there are no boundary surface integrals ds, as we require
# the normal of the variational derivative and test function to vanish 
# at the boundary.


# Define trial and test functions on the mixed space

(u, rho) = TrialFunctions(W)
(dFdu_vec, dFdrho) = TestFunctions(W)

#Define varitional derivatives

(u0,rho0) = split(w0)

dHdu0 = u0/r_0
dHdrho0 = rho0/r_0


# Define discrete divergence operator
def div_u(u, p):
	return (dot(u, grad(p)))*dx(domain=p.ufl_domain()) + (jump(p)*dot((u('-')*(1-theta)+u('+')*theta), n('-')))*dS


L0 = (dot(u0, dFdu_vec) + rho0*dFdrho )*dx
L1 = -div_u(dFdu_vec, r_0 * dHdrho0)
L2 = div_u(dHdu0, r_0 * dFdrho)
#L3 = ( dr_0 * dFdrho * dHdu0[0] - dr_0 * dHdrho0 * dFdu_vec[0] )*dx(degree = quad_degree)

L = L0 + 0.5 * dt * ( L1 + L2   )

a = derivative(L0 - 0.5 * dt * ( L1  + L2 ), w0)

# Storage for visualisation
outfile = File('./Results/compressible_stratified_results.pvd')

u0,rho0 = w0.split()

u0.rename("Velocity")
rho0.rename("Density")


# Output initial conditions
outfile.write(u0,rho0, time = t)


out = Function(W)
# File for energy output
E_file = open('./Results/energy.txt', 'w')






# Solve loop

while (t < end):
    # Update time
    t+= dt
 
    solve(a == L, out, solver_parameters={'ksp_rtol': 1e-14})
    u, rho = out.split()

    # Assign appropriate name in results file
    u.rename("Velocity")
    rho.rename("Density")
 
    # Output results
    #outfile.write(u, rho, time =t)
 
    # Assign output as previous timestep for next time update
    u0.assign(u)
    rho0.assign(rho)

    # Assemble initial energy
    E = assemble ( (0.5/r_0*(inner(u0,u0) + rho0**2))*dx(degree = quad_degree))
 
    E_file.write('%-10s %-10s\n' % (t,abs((E-E0)/E0)))
    # Print time and energy drift, drift should be around machine precision.
    print "At time %g, energy drift is %g" % (t, E-E0)


# Create analytic solutions for error analysis
exact_rho= Function(R)
exact_rho.interpolate(Expression("exp(-0.5* N_sq * x[0])*((N_sq / (2 * sigma))* sin ( m * x[0])+ m / sigma * cos ( m * x[0]))*cos(sigma*(t+0.125))", N_sq = N_sq, m = m , sigma = sigma , t = t))


# Print error for density
error_rho = errornorm(rho, exact_rho,  norm_type='L2')
print "At time %g, l2 error in density is %g" % (t, error_rho)

exact_u = Function(V)
exact_u.interpolate(Expression("exp(-0.5* N_sq * x[0])* sin ( m * x[0])*sin(sigma*(t+0.125))", N_sq = N_sq, m = m , sigma = sigma , t = t))


# Print error for velocity
error_u = errornorm(u, exact_u,  norm_type='L2')
print "At time %g, l2 error in x-velocity is %g" % (t, error_u)

#localenergyfile = File('./Results/local_energy.pvd')

#energy_local = Function(R)
#energy_local.interpolate(dot(u,u)+rho*rho)

#localenergyfile.write(energy_local, time=t)

# Close energy write
E_file.close()
