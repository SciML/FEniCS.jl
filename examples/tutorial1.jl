module Tutorial1
using FEniCS

mesh = UnitSquareMesh(8,8)
V = FunctionSpace(mesh,"P",1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]", degree=2)
u = TrialFunction(V)
bc1 = DirichletBC(V,u_D, "on_boundary")
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u),grad(v))*dx
L = f*v*dx
U=FeFunction(V)
lvsolve(a,L,U,bc1) #linear variational solver
#lvsolve(a,L,U,bc1,solver_parameters=Dict("linear_solver"=>"lu"),form_compiler_parameters=Dict("optimize"=>true))
error_L2= errornorm(u_D, U, norm="L2")

#get_array(L) #this returns an array for the stiffness matrix
#get_array(U) #this returns an array for the solution values
#vtkfile = File("poisson/solution.pvd")
#vtkfile << U.pyobject #exports the solution to a vtkfile
File("poisson/solutionnew.pvd",U)
vertex_values_u_D = compute_vertex_values(u_D,mesh)
vertex_values_u = compute_vertex_values(U,mesh)
error_max = maximum(abs.(vertex_values_u_D-vertex_values_u))
end#module
