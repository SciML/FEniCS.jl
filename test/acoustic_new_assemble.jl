using FEniCS
using PyCall

@pyimport fenics

c = 5000
#problem variables
dt = 0.000004;
t = 0;
T = 0.004;

mesh = RectangleMesh(Point([-2., -2.]),Point([2., 2.]),80,80)
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
#bcs = [BCL,BCD,BCT]

#bc = DirichletBC(V, 0, "on_boundary")
#A, b = assemble_system(a, L, bcs)

u=FEniCS.Function(V)

while t <= T
    delta.pyobject[:t] = t
    #[apply(bcc,b) for bcc in bcs]
    lvsolve(a,L,u,bcs) #linear variational solver
    assign(u0,u1)
    assign(u1,u)
    t +=dt
    #fenics.plot(u.pyobject,title="Acoustic wave Equation")#,mode="auto")

end
true
