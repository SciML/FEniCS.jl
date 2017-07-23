
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

@fenicsclass Expression
Expression(cppcode::String;element=nothing, cell=nothing, domain=nothing, degree=nothing, name=nothing, label=nothing, mpi_comm=nothing) = Expression(fenics.Expression(cppcode=cppcode,
element=element,cell=cell, domain=domain, degree=degree, name=name, label=label, mpi_comm=mpi_comm))
export Expression



#do these all need to be Args or Expression??
inner(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.inner(u.pyobject, v.pyobject))
grad(u::Union{Expression,Argument}) = Expression(fenics.grad(u.pyobject))
nabla_grad(u::Argument) = Expression(fenics.nabla_grad(u.pyobject))
#TODO :need to add dot , etc
export inner,grad, nabla_grad


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
*(expr::Union{Expression,Argument}, expr2::Union{Expression,Argument}) = Expression(expr.pyobject[:__mul__](expr2.pyobject) )
+(measure1::Measure, measure2::Measure) = Form(measure1.pyobject[:__add__](measure2.pyobject) ) #does this need to be expr?
#do we need to export + *?


@fenicsclass Matrix
assemble(assembly_item::Form)=Matrix(fenics.assemble(assembly_item.pyobject)) #this gives as PETScMatrix of appopriate dimensions
export assemble

@fenicsclass BoundaryCondition
DirichletBC(V::FunctionSpace,g,sub_domain,method="topological",check_midpoint=true)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g,sub_domain))#,method=method,check_midpoint=check_midpoint))#look this up with example
#sub_domain seems to have an error I commented kwargs out
export DirichletBC
@fenicsclass sub_domain
#https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/compilemodules/subdomains/CompiledSubDomain.html
CompiledSubDomain(cppcode::String)  = sub_domain(fenics.CompiledSubDomain(cppcode))
export CompiledSubDomain
