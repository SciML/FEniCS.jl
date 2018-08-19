using FEniCS

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)
#in the original tutorial, these are calculated via SymPy
#although it has no actual FEniCS code in calculating those, so I have simply
#stated the result
u_code = "x[0]+2*x[1]+1"
f_code = "-10*x[0]-20*x[1]-10"

calc(var) = 1+var*var

u_D = Expression(u_code, degree=2)
bc = DirichletBC(V, u_D, "on_boundary")
u=FeFunction(V)
v = TestFunction(V)
f = Expression(f_code, degree=2)
F = calc(u)*dot(grad(u), grad(v))*dx - f*v*dx
nlvsolve(F,u,bc)
true
