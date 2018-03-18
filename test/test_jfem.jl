
using FEniCS

mesh = UnitSquareMesh(8, 8)


V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)

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

true
