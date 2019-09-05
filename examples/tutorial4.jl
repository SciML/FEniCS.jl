module Tutorial4
using FEniCS


T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh(Point([-2.0, -2.]), Point([2., 2.]), nx, ny)
V = FunctionSpace(mesh, "P", 1)


bc = DirichletBC(V, Constant(0), "on_boundary")

# Define initial value
u_0 = Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))",
                 degree=2, a=5)

u_n = interpolate(u_0, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Create VTK file for saving solution
#vtkfile = File("heat_gaussian/solution.pvd")

# Time-stepping
u=FeFunction(V)
global t = 0
#fig = PyPlot.figure()
for n = 0:(num_steps-1)
    # Update current time
    global t += dt

    # Compute solution
    lvsolve(a,L,u,bc)


    # Save to file and plot solution
    # vtkfile << (u.pyobject, t)

    # Update previous solution
    assign(u_n,u)
end
end#module
