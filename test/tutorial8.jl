using FEniCS
using PyCall
#using ProgressMeter # for native julia progress

@pyimport fenics


T = 0.05            # final time (increase to 5 for full problem)
num_steps = 50   # number of time steps (#increase to 5000 for full stability)
dt = T / num_steps # time step size
mu = 0.001         # dynamic viscosity
rho = 1            # density

# Create mesh
channel = Rectangle(Point([0.0, 0.0]), Point([2.2, 0.41]))
cylinder = Circle(Point([0.2, 0.2]), 0.05)
domain = channel - cylinder
mesh = generate_mesh(domain, 64)

# Define function spaces
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)

# Define boundaries
inflow   = "near(x[0], 0)"
outflow  = "near(x[0], 2.2)"
walls    = "near(x[1], 0) || near(x[1], 0.41)"
cylinder = "on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3"

# Define inflow profile
inflow_profile = ("4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")

# Define boundary conditions
bcu_inflow = DirichletBC(V, Expression(fenics.Expression(inflow_profile, degree=2)), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = FeFunction(V)
u_  = FeFunction(V)
p_n = FeFunction(Q)
p_ = FeFunction(Q)


U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)


# Define symmetric gradient
function epsilon(u)
    return sym(nabla_grad(u))
end
# Define stress tensor
function sigma(u, p)
    return 2*mu*epsilon(u) - p*Identity(len(u))
end

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx  + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx + inner(sigma(U, p_n), epsilon(v))*dx  + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds   - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)


# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[apply(bc,A1) for bc in bcu]
[apply(bc,A2) for bc in bcp]

# Create XDMF files for visualization output
xdmffile_u = XDMFFile("navier_stokes_cylinder/velocity.xdmf")
xdmffile_p = XDMFFile("navier_stokes_cylinder/pressure.xdmf")

# Create time series (for use in reaction_system.py)
timeseries_u = TimeSeries("navier_stokes_cylinder/velocity_series")
timeseries_p = TimeSeries("navier_stokes_cylinder/pressure_series")

# Save mesh to file (for use in reaction_system.py)
File("navier_stokes_cylinder/cylinder.xml.gz",mesh)




global t = 0
#@showprogress 1 "Solving..." for n=0:num_steps #for native julia progress
for n =0:(num_steps-1)

    # Update current time
    global t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [apply(bc,b1) for bc in bcu]
    solve(A1, vector(u_), b1, "bicgstab", "hypre_amg")

    b2 = assemble(L2)
    [apply(bc,b2) for bc in bcp]
    solve(A2, vector(p_), b2, "bicgstab", "hypre_amg")

    b3 = assemble(L3)
    solve(A3, vector(u_), b3, "cg", "sor")

    write(xdmffile_u,u_, t)
    write(xdmffile_p,p_, t)

    # Save nodal values to file
    store(timeseries_u,vector(u_), t)
    store(timeseries_p,vector(p_), t)
    #update values
    assign(u_n,u_)
    assign(p_n,p_)
end
print("Tutorial 8 finished")

true
