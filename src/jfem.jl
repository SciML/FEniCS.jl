
#These are the commands to define the Fem class, and assemble the Matrix in Julia
#full documentation of the API from FEniCS can be found in the link below
#http://fenics.readthedocs.io/projects/ufl/en/latest/api-doc/ufl.html
#Tests for these can be found in the test_jfem.jl file.


using FEniCS

@fenicsclass FunctionSpace

FunctionSpace(mesh::Mesh, family::Union{String,Symbol}, degree::Int) = FunctionSpace(fenics.FunctionSpace(mesh.pyobject, family, degree))
#add more functionspace functions?
export FunctionSpace

@fenicsclass Argument
Argument(V,number,part::Union{String,Symbol} = nothing) = Argument(fenics.Argument(V.pyobject, number, part=part))
TrialFunction(V::FunctionSpace) = Argument(fenics.TrialFunction(V.pyobject))
TestFunction(V::FunctionSpace) = Argument(fenics.TestFunction(V.pyobject))

export Argument
export TrialFunction
export TestFunction


@fenicsclass Constant

Constant(x::Real) = Constant(fenics.Constant(x, name="Constant($x)"))
export Constant

@fenicsclass Function
Function(V::FunctionSpace) = Function(fenics.Function(V.pyobject))
export Function

@fenicsclass Expression
Expression(cppcode::String;element=nothing, cell=nothing, domain=nothing, degree=nothing, name=nothing, label=nothing, mpi_comm=nothing) = Expression(fenics.Expression(cppcode=cppcode,
element=element,cell=cell, domain=domain, degree=degree, name=name, label=label, mpi_comm=mpi_comm))
#do these all need to be Args or Expression??
inner(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.inner(u.pyobject, v.pyobject))
outer(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.outer(u.pyobject, v.pyobject))
dot(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.dot(u.pyobject, v.pyobject))
grad(u::Union{Expression,Argument}) = Expression(fenics.grad(u.pyobject))
nabla_grad(u::Argument) = Expression(fenics.nabla_grad(u.pyobject))
cross(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.cross(u.pyobject, v.pyobject))
export Expression,inner,grad, nabla_grad,outer,dot,cross
#TODO : CHECK OUTER ,DOT (needs linalg.dot)

@fenicsclass Measure
#cant find the docs for these.
dx = Measure(fenics.dx)
ds = Measure(fenics.ds)
dS = Measure(fenics.dS)
dP = Measure(fenics.dP)
export dx, ds,dS,dP


#https://github.com/FEniCS/ufl/blob/master/ufl/measure.py
@fenicsclass Form
*(expr::Union{Expression,Argument}, measure::Measure) = Form(measure.pyobject[:__rmul__](expr.pyobject) )
*(expr::Union{Expression,Argument,Constant}, expr2::Union{Expression,Argument,Constant}) = Expression(expr.pyobject[:__mul__](expr2.pyobject) )
+(measure1::Measure, measure2::Measure) = Form(measure1.pyobject[:__add__](measure2.pyobject) ) #does this need to be expr?
#do we need to export + *?


@fenicsclass Matrix
#assemble(assembly_item::Form)=Matrix(fenics.assemble(assembly_item.pyobject))
assemble(assembly_item::Form;tensor=nothing, form_compiler_parameters=nothing, add_values=false, finalize_tensor=true, keep_diagonal=false, backend=nothing) = Matrix(
fenics.assemble(assembly_item.pyobject,tensor=tensor,form_compiler_parameters=form_compiler_parameters,add_values=add_values,finalize_tensor=finalize_tensor,keep_diagonal=keep_diagonal,backend=backend))#this gives as PETScMatrix of appopriate dimensions
export assemble

#https://fenicsproject.org/olddocs/dolfin/1.6.0/python/programmers-reference/cpp/fem/DirichletBC.html
@fenicsclass BoundaryCondition
DirichletBC(V::FunctionSpace,g::Expression,sub_domain,method="topological",check_midpoint=true)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g.pyobject,method=method,check_midpoint=check_midpoint))#look this up with example also removed type from g(Should be expression)
export DirichletBC
DirichletBCtest(V::FunctionSpace,g,sub_domain)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g,sub_domain))#look this up with example also removed type from g(Should be expression)

#this class wont be exported until the function boundary has been "fixed"
@fenicsclass SubDomain
#https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/compilemodules/subdomains/CompiledSubDomain.html
CompiledSubDomain(cppcode::String)  = SubDomain(fenics.CompiledSubDomain(cppcode))
export CompiledSubDomain
SubDomain(map_tol::Float64=1e-10)=SubDomain(fenics.SubDomain(map_tol))
export SubDomain
#@fenicsclass MeshFunction
MeshFunction(tp::String;mesh::Mesh=nothing,dim::Int=nothing,filename::String=nothing)=fenics.MeshFunction(tp=tp,mesh=mesh.pyobject,dim=dim,filename=filename)


"""
 For a full list of supported arguments, and their usage
please refer to http://matplotlib.org/api/pyplot_api.html
not all kwargs have been imported. Should you require any that are not imported
open as issue, and I will attempt to add them.
"""
Plot(in_plot::Union{Mesh,FunctionSpace};alpha=1,animated=false,antialiased=true,color="grey"
,dash_capstyle="butt",dash_joinstyle="miter",dashes="",drawstyle="default",fillstyle="full",label="s",linestyle="solid",linewidth=1
,marker="",markeredgecolor="grey",markeredgewidth="",markerfacecolor="grey"
,markerfacecoloralt="grey",markersize=1,markevery="none",visible=true,title="") =fenics.common[:plotting][:plot](in_plot.pyobject,
alpha=alpha,animated=animated,antialiased=antialiased,color=color,dash_capstyle=dash_capstyle,dash_joinstyle=dash_joinstyle
,dashes=dashes,drawstyle=drawstyle,fillstyle=fillstyle,label=label,linestyle=linestyle,linewidth=linewidth,marker=marker,markeredgecolor=markeredgecolor
,markeredgewidth=markeredgewidth,markerfacecolor=markerfacecolor,markerfacecoloralt=markerfacecoloralt,markersize=markersize,markevery=markevery
,visible=visible,title=title)#the first is the keyword argument, the second is the value
export Plot
x123=SubDomain()

@fenicsclass solve
solve(A,x,b)=fenics.solve(A,x,b)


A_15 = [[1,2],[3,4]]
B = [[3,7]]
x = [[]]
#solve(A_15,x,B)
#print(typeof(A_15))
#print(typeof(B))
#print(typeof(x))

mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
#u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

#def boundary(x, on_boundary):
#    return on_boundary

#bc = DirichletBC(V, u_D,6)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx
mf = fenics.MeshFunction("size_t",mesh.pyobject)
# Compute solution
#u = Function(V)
#solve(a == L, u, bc)
#fenics.on_boundary
#DirichletBCtest(V,1,mf)
bc1 = DirichletBCtest(V, 1, "on_boundary")
