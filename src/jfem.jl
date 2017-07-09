
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

TrialFunction(V::FunctionSpace) = Argument(fenics.TrialFunction(V.pyobject))
TestFunction(V::FunctionSpace) = Argument(fenics.TestFunction(V.pyobject))

export TrialFunction
export TestFunction

@fenicsclass Constant

Constant(x::Real) = Constant(fenics.Constant(x, name="Constant($x)"))
export Constant

@fenicsclass Expression
#Expression(expr::String) = Expression(fenics.Expression(expr))

#do these all need to be Args or Expression??
inner(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.inner(u.pyobject, v.pyobject))
grad(u::Union{Expression,Argument}) = Expression(fenics.grad(u.pyobject))
nabla_grad(u::Argument) = Expression(fenics.nabla_grad(u.pyobject))
export inner,grad, nabla_gra


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
