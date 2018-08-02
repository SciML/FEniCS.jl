using FEniCS
using PyCall

@pyimport fenics


L = 1; W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta^2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = BoxMesh(fenics.Point(0, 0, 0), fenics.Point(L, W, W), 10, 3, 3)
V = VectorFunctionSpace(mesh, "P", 1)
c = Constant((0,0,0))
bc = DirichletBC(V, c, "on_boundary && x[1]<1E-14")
tol = 1E-14

#bc1 = fenics.DirichletBC(V.pyobject, c.pyobject, "on_boundary and x[1]<tol",tol=1E-14)
function epsilon(u)
     0.5*(nabla_grad(u)+Transpose(nabla_grad(u)))
end

function sigma(u)
     lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)
end

u = TrialFunction(V)
d = geometric_dimension(u)# space dimension
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds
u=FeFunction(V)
lvsolve(a,L,u,bc)
s = sigma(u) - (1.0/3)*tr(sigma(u))*Identity(d)
#fenics.plot(U.pyobject,mode="displacement")
von_Mises = sqrt(3.0/2*inner(s, s))
V = FunctionSpace(mesh, "P", 1)

von_Mises = project(von_Mises, V)
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
print("Tutorial 6 finished")
true
