
using FEniCS

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)

A_1 = inner(grad(u),grad(v))*dx
B = assemble(A_1)
ng = nabla_grad(u)

measure_x = dx
measure_s = ds
measure_S = dS
measure_P = dP

meas_add1 = dx+ds
meas_add2 = dS+dP

true
