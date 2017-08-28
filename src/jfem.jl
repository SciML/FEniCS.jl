#These are the commands to define the Fem class, and assemble the Matrix in Julia
#full documentation of the API from FEniCS can be found in the link below
#http://fenics.readthedocs.io/projects/ufl/en/latest/api-doc/ufl.html
#Tests for these can be found in the test_jfem.jl file.


using FEniCS

#https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/cpp/fem/FiniteElement.html?highlight=finiteelement
@fenicsclass FiniteElement
FiniteElement(family::Union{String,Symbol},cell=nothing,degree=nothing;form_degree=nothing,quad_scheme=nothing)=FiniteElement(fenics.FiniteElement(family,cell,degree,form_degree=form_degree,
quad_scheme=quad_scheme))
export FiniteElement


@fenicsclass FunctionSpace

FunctionSpace(mesh::Mesh, family::Union{String,Symbol}, degree::Int) = FunctionSpace(fenics.FunctionSpace(mesh.pyobject, family, degree))
FunctionSpace(mesh::Mesh,element::FiniteElement)=FunctionSpace(fenics.FunctionSpace(mesh.pyobject,element.pyobject))
export FunctionSpace

@fenicsclass Argument
Argument(V,number,part::Union{String,Symbol} = nothing) = Argument(fenics.Argument(V.pyobject, number, part=part))
TrialFunction(V::FunctionSpace) = Argument(fenics.TrialFunction(V.pyobject))
TestFunction(V::FunctionSpace) = Argument(fenics.TestFunction(V.pyobject))

export Argument,TrialFunction,TestFunction

@fenicsclass Constant
Constant(x::Real) = Constant(fenics.Constant(x, name="Constant($x)"))
export Constant

@fenicsclass Function
Function(V::FunctionSpace) = Function(fenics.Function(V.pyobject))
export Function

#assigns the computed value
assign(object::Function, solution::Function) = fenicspycall(object, :assign, solution.pyobject)
export assign

@fenicsclass Expression
Expression(cppcode::String;element=nothing, cell=nothing, domain=nothing, degree=nothing, name=nothing, label=nothing, mpi_comm=nothing) = Expression(fenics.Expression(cppcode=cppcode,
element=element,cell=cell, domain=domain, degree=degree, name=name, label=label, mpi_comm=mpi_comm))

inner(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.inner(u.pyobject, v.pyobject))
outer(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.outer(u.pyobject, v.pyobject))
dot(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.dot(u.pyobject, v.pyobject))
grad(u::Union{Expression,Argument}) = Expression(fenics.grad(u.pyobject))
nabla_grad(u::Argument) = Expression(fenics.nabla_grad(u.pyobject))
cross(u::Union{Expression,Argument}, v::Union{Expression,Argument}) = Expression(fenics.cross(u.pyobject, v.pyobject))
export Expression,inner,grad, nabla_grad,outer,dot,cross

@fenicsclass Measure
dx = Measure(fenics.dx)
ds = Measure(fenics.ds)
dS = Measure(fenics.dS)
dP = Measure(fenics.dP)
export dx, ds,dS,dP


#https://github.com/FEniCS/ufl/blob/master/ufl/measure.py
@fenicsclass Form
*(expr::Union{Expression,Argument}, measure::Measure) = Form(measure.pyobject[:__rmul__](expr.pyobject) )
*(expr::Union{Expression,Argument,Constant,Form}, expr2::Union{Expression,Argument,Constant,Form}) = Expression(expr.pyobject[:__mul__](expr2.pyobject) )
*(expr::Float64, expr2::Union{Expression,Argument,Constant}) = Expression(expr2.pyobject[:__mul__](expr) )
*(expr::Union{Expression,Argument,Constant}, expr2::Float64) = Expression(expr.pyobject[:__mul__](expr2) )
+(expr::Union{Expression,Argument,Constant,Measure,Form,Function}, expr2::Union{Expression,Argument,Constant,Measure,Form,Function}) = Expression(expr.pyobject[:__add__](expr2.pyobject) )
-(expr::Union{Expression,Argument,Constant,Measure,Form}, expr2::Union{Expression,Argument,Constant,Measure,Form}) = Expression(expr.pyobject[:__sub__](expr2.pyobject) )

#return form for LHS
lhs(L)=Form(fenics.lhs(L.pyobject))
#return form for RHS
rhs(L)=Form(fenics.rhs(L.pyobject))

export lhs,rhs

#this assembles the matrix from a fenics form
@fenicsclass Matrix
assemble(assembly_item::Union{Form,Expression};tensor=nothing, form_compiler_parameters=nothing, add_values=false, finalize_tensor=true, keep_diagonal=false, backend=nothing) = Matrix(fenics.assemble(assembly_item.pyobject,
tensor=tensor,form_compiler_parameters=form_compiler_parameters,add_values=add_values,finalize_tensor=finalize_tensor,keep_diagonal=keep_diagonal,backend=backend))
export assemble
#I have changed this to Function+Form

#https://fenicsproject.org/olddocs/dolfin/1.6.0/python/programmers-reference/cpp/fem/DirichletBC.html
@fenicsclass sub_domain
@fenicsclass BoundaryCondition

"""
DirichletBC works in a slightly altered way to the original FEniCS.
Instead of defining the function with the required return command,
we define the return command directly (
ie instead of function(x)
                 return x
we simple write "x"
"""
DirichletBC(V::FunctionSpace,g,sub_domain)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g.pyobject,sub_domain))#look this up with example also removed type from g(Should be expression)
export DirichletBC
#https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/compilemodules/subdomains/CompiledSubDomain.html
CompiledSubDomain(cppcode::String)  = sub_domain(fenics.CompiledSubDomain(cppcode))
export CompiledSubDomain
#This function provides a SubDomain which picks out the boundary of a mesh, and provides a convenient way to specify boundary conditions on the entire boundary of a mesh.
#so it can be used in the place of "on_boundary" etc.
DomainBoundary()=fenics.DomainBoundary()
export DomainBoundary

"""
 For a full list of supported arguments, and their usage
please refer to http://matplotlib.org/api/pyplot_api.html
not all kwargs have been imported. Should you require any that are not imported
open as issue, and I will attempt to add them.
"""
Plot(in_plot::Union{Mesh,FunctionSpace,Function,Geometry};alpha=1,animated=false,antialiased=true,color="grey"
,dash_capstyle="butt",dash_joinstyle="miter",dashes="",drawstyle="default",fillstyle="full",label="s",linestyle="solid",linewidth=1
,marker="",markeredgecolor="grey",markeredgewidth="",markerfacecolor="grey"
,markerfacecoloralt="grey",markersize=1,markevery="none",visible=true,title="") =fenics.common[:plotting][:plot](in_plot.pyobject,
alpha=alpha,animated=animated,antialiased=antialiased,color=color,dash_capstyle=dash_capstyle,dash_joinstyle=dash_joinstyle
,dashes=dashes,drawstyle=drawstyle,fillstyle=fillstyle,label=label,linestyle=linestyle,linewidth=linewidth,marker=marker,markeredgecolor=markeredgecolor
,markeredgewidth=markeredgewidth,markerfacecolor=markerfacecolor,markerfacecoloralt=markerfacecoloralt,markersize=markersize,markevery=markevery
,visible=visible,title=title)#the first is the keyword argument, the second is the value
#export Plot
