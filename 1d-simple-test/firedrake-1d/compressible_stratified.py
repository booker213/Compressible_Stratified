from firedrake import *

try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

# Script to solve compressible stratified waves in upto three dimensions
# c^2_0 and gravity will be considered to be scaled to be equal to 1.

# Test problem will be the one-dimensional waves considered in the MATLAB scripts.

# *_n1 represents variables at time n + 1
# *_n represents variables at time n 

# Create mesh
# Current mesh is a unit line  with Nx elements.

Nx = 4
mesh = UnitIntervalMesh(Nx)

# Declare timestep
dt = 1./( pow(Nx,2))

# Declare initial and end time
# Period of waves considered in test is 1s
# We will consider 3 periods initially
t = 0.
end =  10.

# Declare order of the basis in the elements
order_basis = 6

# The projection of the variational derivatives requires a higher degree
# of quadrature than the automatic procedure provides. Here we define 
# a sufficiently high quadrative degree.

quad_degree = 2 * order_basis

# Declare flux indicator function
theta = Constant(0.6)

# Declare constants 
N_sq = Constant(3.0) # Buoyancy Frequency
m = 2 * pi  # Wave number
sigma = pow( 0.25 * pow(N_sq, 2) + pow(m, 2),0.5) # Wave Frequency


# Declare function spaces on the mesh
# and create mixed function space
# for coupled Poisson Bracket.

V = VectorFunctionSpace(mesh, "DG", order_basis)
R = FunctionSpace(mesh, "DG", order_basis)

W = V*R*V*R




# Define Background Density
rho_0_project = Function(R)
rho_0_project.interpolate(Expression( "exp(-3.0*x[0])" ))






# Initial conditions
# Here time n = 0

# Function space
w_n = Function(W)


# Interpolate expressions
# Use .split() to assign values to the components of the mixed vector
u_n,rho_n, dHdu_n, dHdrho_n = w_n.split()

u_n.interpolate(Expression("exp(-0.5* N_sq * x[0])* sin ( m * x[0])*sin(sigma*0.125)", N_sq = N_sq, m = m , sigma = sigma ))

rho_n.interpolate(Expression("exp(-0.5* N_sq * x[0])*((N_sq / (2 * sigma))* sin ( m * x[0])+ m / sigma * cos ( m * x[0]))*cos(sigma*0.125)", N_sq = N_sq, m = m , sigma = sigma ))

# Project initial variational derivatives

q = TrialFunction(V)
q_test = TestFunction(V)

a_u_project = dot(q ,  q_test) * dx(degree = quad_degree)
L_u_project = dot (u_n / rho_0_project ,  q_test) * dx(degree = quad_degree)

b = TrialFunction(R)
b_test = TestFunction(R)

a_rho_project = b * b_test * dx(degree = quad_degree)
L_rho_project = rho_n / rho_0_project * b_test * dx(degree = quad_degree)


solve ( a_u_project == L_u_project , dHdu_n , solver_parameters={'ksp_rtol': 1e-14} )
solve ( a_rho_project == L_rho_project , dHdrho_n, solver_parameters={'ksp_rtol': 1e-14} )



# Assemble initial energy
E0 = assemble ( (0.5/rho_0_project*(inner(u_n,u_n) + rho_n**2))*dx(degree = quad_degree))

print "The initial energy is %g " % (E0)

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

(u_n1, rho_n1, dHdu_n1, dHdrho_n1) = TrialFunctions(W)
(dFdu_vec, dFdrho, dFdu_project, dFdrho_project) = TestFunctions(W)

#  Use split(*) to use initial conditions in the variational form

(u_n,rho_n, dHdu_n, dHdrho_n) = split(w_n)


# Define discrete divergence operator
def div_u(u, p):
	return (dot(u, grad(p)))*dx(domain=p.ufl_domain()) + (jump(p)*dot((u('-')*(1-theta)+u('+')*theta), n('-')))*dS


L0 = (dot(u_n, dFdu_vec) + rho_n*dFdrho )*dx
L1 = -div_u(dFdu_vec, rho_0_project * dHdrho_n)
L2 = div_u(dHdu_n, rho_0_project * dFdrho)


L = L0 + 0.5 * dt * ( L1 + L2   )

a1 = ( dot(dHdu_n1 , dFdu_project) - dot(u_n1 / rho_0_project ,  dFdu_project) ) * dx(degree = quad_degree)
a2 = ( dHdrho_n1 * dFdrho_project - (rho_n1 / rho_0_project )* dFdrho_project ) *dx(degree = quad_degree)

a = derivative(L0 - 0.5 * dt * ( L1  + L2 ), w_n) + a1 + a2
 


# Storage for visualisation
outfile = File('./Results/compressible_stratified_results.pvd')

# Use .split() to access the data.
u_n, rho_n, dHdu_n, dHdrho_n = w_n.split()

u_n.rename("Velocity")
rho_n.rename("Density")


# Output initial conditions
outfile.write(u_n,rho_n, time = t)
# Assign array storage for matplotlib plotting
#velocity = []
density = []

out = Function(W)
# File for energy output
E_file = open('./Results/energy.txt', 'w')


# Solve loop

while (t < end):
    # Update time
    t+= dt
 
    solve(a == L, out, solver_parameters={'ksp_rtol': 1e-14})
    u_n1, rho_n1, dHdu_n1, dHdrho_n1 = out.split()

    # Assign appropriate name in results file
    u_n1.rename("Velocity")
    rho_n1.rename("Density")
 
    # Output results
    outfile.write(u_n1, rho_n1, time =t)
    
    # Write  plotting data to array
    #velocity.append(Function(u_n1)) # Currently not supported for vector objects
    density.append(Function(rho_n1))
    
    # Assign output as previous timestep for next time update
    u_n.assign(u_n1)
    rho_n.assign(rho_n1)
    dHdu_n.assign(dHdu_n1)
    dHdrho_n.assign(dHdrho_n1)
 
    # Assemble initial energy
    E = assemble ( (0.5/rho_0_project*(inner(u_n,u_n) + rho_n**2))*dx(degree = quad_degree))
 
    E_file.write('%-10s %-10s\n' % (t,abs((E-E0)/E0)))
    # Print time and energy drift, drift should be around machine precision.
    print "At time %g, energy drift is %g" % (t, E-E0)




# Use matplotlib to plot scalar results
try:
 plot(density)
except Exception as e:
 warning("Cannot plot figure. Error msg: '%s'" % e.message)
try:
 plt.show()
except Exception as e:
 warning("Cannot show figure. Error msg: '%s'" % e.message)

# Create analytic solutions for error analysis
exact_rho= Function(R)
exact_rho.interpolate(Expression("exp(-0.5* N_sq * x[0])*((N_sq / (2 * sigma))* sin ( m * x[0])+ m / sigma * cos ( m * x[0]))*cos(sigma*(t+0.125))", N_sq = N_sq, m = m , sigma = sigma , t = t))


# Print error for density
error_rho = errornorm(rho_n1, exact_rho,  norm_type='L2')
print "At time %g, l2 error in density is %g" % (t, error_rho)

exact_u = Function(V)
exact_u.interpolate(Expression("exp(-0.5* N_sq * x[0])* sin ( m * x[0])*sin(sigma*(t+0.125))", N_sq = N_sq, m = m , sigma = sigma , t = t))


# Print error for velocity
error_u = errornorm(u_n1, exact_u,  norm_type='L2')
print "At time %g, l2 error in x-velocity is %g" % (t, error_u)



# Close energy write
E_file.close()
