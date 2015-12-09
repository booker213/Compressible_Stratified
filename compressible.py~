from firedrake import *

#Create Mesh
m=16
mesh = UnitSquareMesh(m, m)
order_basis = 0 
#Declare timestep 
timestep = 1/16
t=0.0
end=1

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

bcv_x = DirichletBC(W.sub(0).sub(0), Constant(0), 1, method="geometric")
bcv_z = DirichletBC(W.sub(0).sub(1), Constant(0), 1, method="geometric")

#outfile = File("./Results/compressible.pvd")

#Define Poisson Bracket
#Possibly as a python definition


#while (t <= end):

