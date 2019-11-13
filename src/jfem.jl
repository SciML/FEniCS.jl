#These are the commands to define the Fem class, and assemble the Matrix in Julia
#full documentation of the API from FEniCS can be found in the link below
#http://fenics.readthedocs.io/projects/UFL/en/latest/api-doc/Expression.html
#Tests for these can be found in the test_jfem.jl file.


using FEniCS

@fenicsclass FunctionSpace

#Use Functionspace for scalar fields
FunctionSpace(mesh::Mesh, family::StringOrSymbol, degree::Int) = FunctionSpace(fenics.FunctionSpace(mesh.pyobject, family, degree))
#Use VectorFunctionSpace for vector fields
VectorFunctionSpace(mesh::Mesh, family::StringOrSymbol,degree::Int) = FunctionSpace(fenics.VectorFunctionSpace(mesh.pyobject, family, degree))

export FunctionSpace, VectorFunctionSpace

@fenicsclass Expression
Argument(V,number,part::StringOrSymbol = nothing) = Expression(fenics.Argument(V.pyobject, number, part=part))
TrialFunction(V::FunctionSpace) = Expression(fenics.TrialFunction(V.pyobject))

TestFunction(V::FunctionSpace) = Expression(fenics.TestFunction(V.pyobject))

function TrialFunctions(V::FunctionSpace)
     vec=fenics.TrialFunctions(V.pyobject)
     expr_vec = [Expression(elem) for elem in vec]
     return expr_vec
end

function TestFunctions(V::FunctionSpace)
     vec=fenics.TestFunctions(V.pyobject)
     expr_vec = [Expression(elem) for elem in vec]
     return expr_vec
end



#Below are attributes for the argument *class*



export Argument,TrialFunction,TrialFunctions,TestFunction,TestFunctions

@fenicsclass Constant
Constant(x::Union{Real,Tuple}) = Expression(fenics.Constant(x, name="Constant($x)"))
export Constant

@fenicsclass FeFunction
function FeFunction(V::FunctionSpace; name::String="")
    if name == ""
        return FeFunction(fenics.Function(V.pyobject))
    else
        return FeFunction(fenics.Function(V.pyobject, name=name))
    end
end
assign(solution1::FeFunction,solution2 )=fenicspycall(solution1,:assign,solution2.pyobject)

function assign(solution::FeFunction, data::AbstractArray)
    solution.pyobject.vector().set_local(data)
end

geometric_dimension(expr::Union{FeFunction,Expression}) = fenicspycall(expr, :geometric_dimension)
export geometric_dimension

function split(fun::FeFunction)
     vec = fenics.split(fun.pyobject)
     expr_vec = [FeFunction(spl) for spl in vec]
     return expr_vec
end

function py_split(fun::FeFunction)
     vec = fun.pyobject.split()
     expr_vec = [FeFunction(spl) for spl in vec]
     return expr_vec
end

export FeFunction,assign,split,py_split

function Expression(cppcode; kw...)
    Expression(fenics.Expression(cppcode; kw...))
end
Identity(dim::Int) = Expression(fenics.Identity(dim))
inner(u::Union{Expression,FeFunction}, v::Union{Expression,FeFunction}) = Expression(fenics.inner(u.pyobject, v.pyobject))
outer(u::Union{Expression,FeFunction}, v::Union{Expression,FeFunction}) = Expression(fenics.outer(u.pyobject, v.pyobject))
dot(u::Union{Expression,FeFunction}, v::Union{Expression,FeFunction}) = Expression(fenics.dot(u.pyobject, v.pyobject))
grad(u::Union{Expression,FeFunction}) = Expression(fenics.grad(u.pyobject))
∇(u::Union{Expression,FeFunction}) = Expression(fenics.grad(u.pyobject))
nabla_grad(u::Union{Expression,FeFunction}) = Expression(ufl.nabla_grad(u.pyobject))
nabla_div(u::Union{Expression,FeFunction}) = Expression(ufl.nabla_div(u.pyobject))
div(u::Union{Expression,FeFunction}) = Expression(fenics.div(u.pyobject))
cross(u::Union{Expression,FeFunction}, v::Union{Expression,FeFunction}) = Expression(fenics.cross(u.pyobject, v.pyobject))
tr(u::Union{Expression,FeFunction}) = Expression(fenics.tr(u.pyobject))
sqrt(u::Union{Expression,FeFunction}) = Expression(fenics.sqrt(u.pyobject))
sym(u::Union{Expression,FeFunction}) = Expression(fenics.sym(u.pyobject))
len(U::Union{Expression,FeFunction}) = length(U.pyobject)

interpolate(solution1::FeFunction,solution2::Expression) = FeFunction(fenicspycall(solution1,:interpolate,solution2.pyobject))

Expression(x::FEniCS.Expression) = convert(Expression,x)
export Expression,Identity,inner,grad, nabla_grad, nabla_div,div, outer,dot,cross, tr, sqrt, sym, len, interpolate
export ∇

#Below are attributes for the Expression and FeFunction types
#Computes values at vertex of a mesh for a given Expression / FeFunction, and returns an array
compute_vertex_values(expr::Expression, mesh::Mesh) = fenicspycall(expr, :compute_vertex_values, mesh.pyobject)
compute_vertex_values(expr::FeFunction, mesh::Mesh) = fenicspycall(expr, :compute_vertex_values, mesh.pyobject)

export compute_vertex_values

@fenicsclass Measure
directional_derivative(solution1::FeFunction,direction)=FeFunction(fenicspycall(solution1,:dx,direction))

# some of these constant are initialized in __init__
export dx,ds,dS,dP,directional_derivative

#https://github.com/FEniCS/Expression/blob/master/Expression/measure.py
@fenicsclass Form

*(expr::Union{Expression,FeFunction}, measure::Measure) = Expression(measure.pyobject.__rmul__(expr.pyobject) )

*(expr::Union{Expression,FeFunction}, expr2::Union{Expression,FeFunction}) = Expression(expr.pyobject.__mul__(expr2.pyobject) )
*(expr::Real, expr2::Union{Expression,FeFunction}) = Expression(expr2.pyobject.__mul__(expr) )
*(expr::Union{Expression,FeFunction}, expr2::Real) = Expression(expr.pyobject.__mul__(expr2) )

+(expr::Union{Expression,FeFunction}, expr2::Real) = Expression(expr.pyobject.__add__(expr2) )
+(expr::Real, expr2::Union{Expression,FeFunction}) = Expression(expr2.pyobject.__add__(expr) )
+(expr::Union{Expression,FeFunction}, expr2::Union{Expression,FeFunction}) = Expression(expr.pyobject.__add__(expr2.pyobject) )

-(expr::Union{Expression,FeFunction}, expr2::Real) = Expression(expr.pyobject.__sub__(expr2) )
-(expr::Real, expr2::Union{Expression,FeFunction}) = -1*(Expression(expr2.pyobject.__sub__(expr) ))
-(expr::Union{Expression,FeFunction}, expr2::Union{Expression,FeFunction}) = Expression(expr.pyobject.__sub__(expr2.pyobject) )

/(expr::Union{Expression,FeFunction}, expr2::Real) = Expression(expr.pyobject.__div__(expr2) )
/(expr::Union{Expression,FeFunction}, expr2::Union{Expression,FeFunction}) = Expression(expr.pyobject.__div__(expr2.pyobject) )

function /(expr::Real,expr2::Union{Expression,FeFunction})
    x = expr2*expr2
    y = x/expr
    z = expr2/y
    return Expression(z)
end

function Transpose(object::Expression)
    x = object.pyobject.T
    y = Expression(x)
    return y
end
export Transpose

"""
rhs(equation::Expression)
Given a combined bilinear and linear form,
    extract the right hand side (negated linear form part).
    Example::
           a = u*v*dx + f*v*dx
           L = rhs(a) -> -f*v*dx
"""
rhs(equation::Expression)=Expression(fenics.rhs(equation.pyobject))
"""
lhs(equation::Expression)
    Given a combined bilinear and linear form,
    extract the left hand side (bilinear form part).
    Example::
        a = u*v*dx + f*v*dx
        a = lhs(a) -> u*v*dx
"""
lhs(equation::Expression)=Expression(fenics.lhs(equation.pyobject))

export lhs, rhs
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
DirichletBC(V::FunctionSpace,g::Number,sub_domain)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g,sub_domain))#look this up with example also removed type from g(Should be expression)
DirichletBC(V::FunctionSpace,g::Tuple,sub_domain)=BoundaryCondition(fenics.DirichletBC(V.pyobject,g,sub_domain))#look this up with example also removed type from g(Should be expression)

export DirichletBC
#https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/compilemodules/subdomains/CompiledSubDomain.html
CompiledSubDomain(cppcode::String)  = sub_domain(fenics.CompiledSubDomain(cppcode))
export CompiledSubDomain

apply(bcs::BoundaryCondition, matrix::Matrix) = fenicspycall(bcs, :apply, matrix.pyobject)
apply(bcs, matrix::Matrix) = fenicspycall(BoundaryCondition(bcs), :apply, matrix.pyobject)

export apply

function assemble_system_julia(a::Expression,L::Expression)
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject)
    A  = A_fenics[:array]()
    b = b_fenics[:array]()
    return A,b
end

function assemble_system_julia(a::Expression,L::Expression,bc)
    bcs_py = [bcs.pyobject for bcs in bc]
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject,bcs_py)
    A  = A_fenics[:array]()
    b = b_fenics[:array]()
    return A,b
end

function assemble_system_julia(a::Expression,L::Expression,bc::BoundaryCondition)
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject,bc.pyobject)
    A  = A_fenics[:array]()
    b = b_fenics[:array]()
    return A,b
end

##add assemble_system to FEniCS objects for solving

function assemble_system(a::Expression,L::Expression)
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject)
    A  = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A,b
end

function assemble_system(a::Expression,L::Expression,bc)
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject,bc)
    A  = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A,b
end

function assemble_system(a::Expression,L::Expression,bc::BoundaryCondition)
    A_fenics,b_fenics = fenics.assemble_system(a.pyobject,L.pyobject,bc)
    A  = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A,b
end


export assemble_system, assemble_system_julia


"""
 For a full list of supported arguments, and their usage
please refer to http://matplotlib.org/api/pyplot_api.html
not all kwargs have been imported. Should you require any that are not imported
open as issue, and I will attempt to add them.
Deprecate this in a future version
"""
Plot(in_plot::Union{Mesh,FunctionSpace,FeFunction};alpha=1,animated=false,antialiased=true,color="grey"
,dash_capstyle="butt",dash_joinstyle="miter",dashes="",drawstyle="default",fillstyle="full",label="s",linestyle="solid",linewidth=1
,marker="",markeredgecolor="grey",markeredgewidth="",markerfacecolor="grey"
,markerfacecoloralt="grey",markersize=1,markevery="none",visible=true,title="") =fenics.common.plotting.plot(in_plot.pyobject,
alpha=alpha,animated=animated,antialiased=antialiased,color=color,dash_capstyle=dash_capstyle,dash_joinstyle=dash_joinstyle
,dashes=dashes,drawstyle=drawstyle,fillstyle=fillstyle,label=label,linestyle=linestyle,linewidth=linewidth,marker=marker,markeredgecolor=markeredgecolor
,markeredgewidth=markeredgewidth,markerfacecolor=markerfacecolor,markerfacecoloralt=markerfacecoloralt,markersize=markersize,markevery=markevery
,visible=visible,title=title)#the first is the keyword argument, the second is the value
#export Plot



@fenicsclass FiniteElement
"""
This is the FiniteElement class. Attributes/function inputs can be found below,
or via the FEniCS documentation.
This is only a basic implementation. Many of the expected methods are still missing
Arguments
|          family (string/symbol)
|             The finite element family
|          cell
|             The geometric cell
|          degree (int)
|             The polynomial degree (optional)
|          form_degree (int)
|             The form degree (FEEC notation, used when field is
|             viewed as k-form)
|          quad_scheme
|             The quadrature scheme (optional)
|          variant
|             Hint for the local basis function variant (optional)
"""
FiniteElement(family::StringOrSymbol,cell=nothing,degree=nothing,form_degree=nothing,quad_scheme=nothing,variant=nothing) =
FiniteElement(fenics.FiniteElement(family=family, cell=cell, degree=degree, form_degree=form_degree, quad_scheme=quad_scheme,variant=variant))
export FiniteElement

"""
Methods for the FiniteElement class
mapping(self)
 |
 |  reconstruct(self, family=None, cell=None, degree=None)
 |      Construct a new FiniteElement object with some properties
 |      replaced with new values.
 |
 |  shortstr(self)
 |      Format as string for pretty printing.
 |
 |  sobolev_space(self)
 |      Return the underlying Sobolev space.
 |
 |  variant(self)
 |
"""

# some of these constant are initialized in __init__
export hexahedron, tetrahedron, quadrilateral, triangle

family(finiteelement::FiniteElement) = fenicspycall(finiteelement, :family)

cell(finiteelement::FiniteElement) = fenicspycall(finiteelement, :cell)

degree(finiteelement::FiniteElement) = fenicspycall(finiteelement, :degree)

#form_degree(finiteelement::FiniteElement) = fenicspycall(finiteelement, :form_degree)

#quad_scheme(finiteelement::FiniteElement) = fenicspycall(finiteelement, :quad_scheme)

variant(finiteelement::FiniteElement) = fenicspycall(finiteelement, :variant)

reconstruct(finiteelement::FiniteElement, ;family=nothing,cell=nothing,degree=nothing) =
FiniteElement(fenicspycall(finiteelement, :reconstruct, family,cell,degree))

sobolev_space(finiteelement::FiniteElement) = fenicspycall(finiteelement, :sobolev_space)
export family, cell, degree, variant, reconstruct, sobolev_space

FacetNormal(mesh::Mesh)  = Expression(fenics.FacetNormal(mesh.pyobject))

export FacetNormal

@fenicsclass MixedElement

function MixedElement(vec::Array{FEniCS.FiniteElementImpl,1})
    pyvec = [elem.pyobject for elem in vec]
    me = MixedElement(fenics.MixedElement(pyvec))
end

export MixedElement

FunctionSpace(mesh::Mesh, element::Union{FiniteElement,MixedElement}) = FunctionSpace(fenics.FunctionSpace(mesh.pyobject, element.pyobject))

export FunctionSpace
