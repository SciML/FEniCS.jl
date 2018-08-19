using FEniCS
using PyCall

@pyimport fenics
c = 5000
#problem variables, have been scaled down for faster test solution
dt = 0.00004;
global t = 0;
T = 0.0004;

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

#bcs = [BCL.pyobject,BCD.pyobject,BCT.pyobject,BCR.pyobject]
bcs = [BCL.pyobject,BCD.pyobject,BCT.pyobject]

#bc = DirichletBC(V, 0, "on_boundary")

u=FeFunction(V)

while t <= T
    A, b = assemble_system(a, L, bcs)
    delta.pyobject[:t] = t
    [apply(bcc,b) for bcc in bcs]
    solve(A, vector(u), b)
    assign(u0,u1)
    assign(u1,u)
    global t +=dt
    #fenics.plot(u.pyobject,title="Acoustic wave Equation")#,mode="auto")

end
println("Acoustic problem finished")
true
