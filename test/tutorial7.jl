using FEniCS
using PyCall

@pyimport fenics


T = 1.0           # final time (increase to 10 for full problem)
num_steps = 50    # number of time steps (increase for stability)
dt = T / num_steps # time step size
mu = 1             # kinematic viscosity
rho = 1            # density

# Create mesh and define function spaces
mesh = UnitSquareMesh(16, 16)
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)

inflow  = "near(x[0], 0)"
outflow = "near(x[0], 1)"
walls = "near(x[1], 0) || near(x[1], 1)"
bcu_noslip  = DirichletBC(V, Constant((0, 0)), walls)
bcp_inflow  = DirichletBC(Q, Constant(8), inflow)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)

bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]


u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

u_n = FeFunction(V)
u_  = FeFunction(V)
p_n = FeFunction(Q)
p_ = FeFunction(Q)

U   = 0.5*(u_n + u)
n   = FacetNormal(mesh)
f   = Constant((0, 0))
k   = Constant(dt)
mu  = Constant(mu)
rho = Constant(rho)

function epsilon(u)
    return sym(nabla_grad(u))
end

function sigma(u, p)
    return 2*mu*epsilon(u) - p*Identity(len(u))
end

F1 = rho*dot((u - u_n) / k, v)*dx + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx  + inner(sigma(U, p_n), epsilon(v))*dx  + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx


A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

[apply(bc,A1) for bc in bcu]
[apply(bc,A2) for bc in bcp]

# Time-stepping
global t = 0
for n = 0:(num_steps-1)

    # Update current time
    global t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [apply(bc,b1) for bc in bcu]

    solve(A1,vector(u_),b1)

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [apply(bc,b2) for bc in bcp]
    solve(A2, vector(p_), b2)

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, vector(u_), b3)

    # Plot solution
    #fenics.plot(u_.pyobject)

    # Compute error
    u_e = Expression(fenics.Expression(("4*x[1]*(1.0 - x[1])", "0"), degree=2))
    u_e = interpolate(u_e, V)

    # Update previous solution
    assign(u_n,u_)
    assign(p_n,p_)
    #println(n)
# Hold plot
end
print("Tutorial 7 finished")
true
