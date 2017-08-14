
using FEniCS

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]",degree=2 )
bc1 = DirichletBC(V,u_D, "on_boundary")
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u),grad(v))*dx
L = f*v*dx
U = FEniCS.Function(V)

lvsolve(a,L,U,bc1)
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
true
