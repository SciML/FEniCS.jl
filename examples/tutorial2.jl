module Tutorial2
using FEniCS
domain = Circle(Point([0.0,0.0]),1)
mesh = generate_mesh(domain,16)
V = FunctionSpace(mesh,"P",2)
w_D = Constant(0)

bc = DirichletBC(V,w_D,"on_boundary")
beta =8
R0 = 0.6
p = Expression("4*exp(-pow(beta, 2)*(pow(x[0], 2) + pow(x[1] - R0, 2)))",
degree=1, beta=beta, R0=R0)
w = TrialFunction(V)
v = TestFunction(V)
a = dot(grad(w), grad(v))*dx
L = p*v*dx
w=FeFunction(V)
lvsolve(a,L,w,bc)
p = interpolate(p,V)

vtkfile_w = File("poisson_membrane/deflection.pvd")
vtkfile_w << w.pyobject
vtkfile_p = File("poisson_membrane/load.pvd")
vtkfile_p << p.pyobject
end#module
