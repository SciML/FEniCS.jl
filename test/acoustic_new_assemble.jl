using FEniCS
using PyCall

@pyimport fenics

c = 5000
#problem variables
dt = 0.000004;
global t = 0;
T = 0.004;

mesh = RectangleMesh(Point([-2., -2.]),Point([2., 2.]),40,40)
V = FunctionSpace(mesh,"Lagrange",1)

# Previous and current solution
u1= interpolate(Constant(0.0), V)
u0= interpolate(Constant(0.0), V)

# Variational problem at each time
u = TrialFunction(V)
v = TestFunction(V)

a = u*v*dx + dt*dt*c*c*inner(grad(u), grad(v))*dx
L = 2*u1*v*dx-u0*v*dx

delta = Expression(fenics.Expression("sin(c*10*t)",degree=2,c=c,t=t))

#Define boundary conditions
top ="on_boundary && near(x[1], 2)"
BCT = DirichletBC(V,Constant(0.0),top)

down ="on_boundary && near(x[1], -2)"
BCD = DirichletBC(V,Constant(0.0),down)

left = "on_boundary && near(x[0],-2)"
BCL = DirichletBC(V,delta,left)

right = "on_boundary && near(x[0],2)"
BCR = DirichletBC(V,Constant(0.0),right)

bcs_dir = [BCL,BCD,BCT,BCR]
bcs_neu = [BCL,BCD,BCT]


u=FeFunction(V)

while t <= T
    delta.pyobject[:t] = t
    lvsolve(a,L,u,bcs_dir) #linear variational solver
    assign(u0,u1)
    assign(u1,u)
    global t +=dt
end
println("Acoustic problem finished")
true
