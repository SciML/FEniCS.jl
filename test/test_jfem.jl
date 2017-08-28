
using FEniCS
using PyCall
@pyimport fenics

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]",degree=2 )
bc = DirichletBC(V,u_D, "on_boundary")
element = FiniteElement("Lagrange",ufl_cell(mesh),1)
#V_new = FunctionSpace(mesh,element)
bc_new = DirichletBC(V,u_D, DomainBoundary())

u = TrialFunction(V)
v = TestFunction(V)
f = FEniCS.Constant(-6.0)
a = dot(grad(u),grad(v))*dx
L = f*v*dx
U = FEniCS.Function(V)
U_new = FEniCS.Function(V)
lvsolve(a,L,U,bc)
lvsolve(a,L,U_new,bc_new)
errornorm(u_D,U_new,norm="L2")
errornorm(u_D,U,norm="L2")
file_sol= File("FEniCS.jl/test_jfemsol.pvd")
file_sol << U.pyobject



A_1 = inner(grad(u),grad(v))*dx
A_2 = dot(grad(u),grad(v))*dx
A_3 = cross(grad(u),grad(v))
A_4 = outer(grad(u),grad(v))
ng = nabla_grad(u)

B_1 = assemble(A_1)
B_2 = assemble(A_2)

measure_x = dx
measure_s = ds
measure_S = dS
measure_P = dP

meas_add1 = dx+ds
meas_add2 = dS+dP

#matrix_U = assemble(U)
matrix_L = assemble(L)
vector_U = vector(U)
genericvector_U = Vector(vector_U)

T = 2.0            # final time
num_steps = 10     # number of time steps
dt = T / num_steps # time step size
alpha = 3          # parameter alpha
beta = 1.2         # parameter beta
t=0

u_D = Expression(fenics.Expression("1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t",
                 degree=2, alpha=alpha, beta=beta, t=0))

bc = DirichletBC(V, u_D, "on_boundary")

u_n = interpolate(u_D, V)
u_1 = TrialFunction(V)
v_1 = TestFunction(V)
f_1= Constant(beta - 2 - 2*alpha)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx

A_1, L_1 = lhs(F), rhs(F)
M = FEniCS.Function(V)

#vtkfile = File("FEniCS.jl/test/location.pvd")

for n =1:num_steps
    t += dt
    u_D.pyobject[:t] =t
    lvsolve(A_1 , L_1, M, bc)
#    vtkfile << (M.pyobject, t) #if you wish to save each time_step
    u_e = interpolate(u_D, V)
    assign(u_n,M)
end

true
