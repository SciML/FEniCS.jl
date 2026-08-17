#These are the commands to define the Fem class, and assemble the Matrix in Julia
#full documentation of the API from FEniCS can be found in the link below
#http://fenics.readthedocs.io/projects/UFL/en/latest/api-doc/Expression.html
#Tests for these can be found in the test_jfem.jl file.

"""
    FunctionSpace

Abstract wrapper for a FEniCS finite-element function space.
"""
abstract type FunctionSpace <: fenicsobject end

struct FunctionSpaceImpl <: FunctionSpace
    pyobject::PyObject
end

FunctionSpace(pyobject::PyObject) = FunctionSpaceImpl(pyobject)

"""
    FunctionSpace(mesh::Mesh, family::StringOrSymbol, degree::Int)

Construct a scalar finite-element function space on `mesh`.

# Arguments

- `mesh`: Mesh on which the space is defined.
- `family`: FEniCS finite-element family, such as `"CG"`.
- `degree`: Polynomial degree of the basis.
"""
function FunctionSpace(mesh::Mesh, family::StringOrSymbol, degree::Int)
    return FunctionSpace(fenics.FunctionSpace(mesh.pyobject, family, degree))
end
"""
    VectorFunctionSpace(mesh::Mesh, family::StringOrSymbol, degree::Int)

Construct a vector-valued finite-element function space on `mesh`.

# Arguments

- `mesh`: Mesh on which the space is defined.
- `family`: FEniCS finite-element family.
- `degree`: Polynomial degree of the basis.
"""
function VectorFunctionSpace(mesh::Mesh, family::StringOrSymbol, degree::Int)
    return FunctionSpace(fenics.VectorFunctionSpace(mesh.pyobject, family, degree))
end

export FunctionSpace, VectorFunctionSpace

"""
    Expression

Abstract wrapper for a FEniCS/UFL symbolic expression.
"""
abstract type Expression <: fenicsobject end

struct ExpressionImpl <: Expression
    pyobject::PyObject
end

Expression(pyobject::PyObject) = ExpressionImpl(pyobject)
"""
    Argument(V, number, part = nothing)

Construct a symbolic UFL argument associated with function space `V`.

# Arguments

- `V`: Function space associated with the argument.
- `number`: Argument number used by UFL.
- `part`: Optional component identifier for mixed spaces.
"""
function Argument(V, number, part::Union{StringOrSymbol, Nothing} = nothing)
    return Expression(fenics.Argument(V.pyobject, number, part = part))
end
"""
    TrialFunction(V::FunctionSpace)

Construct a symbolic trial function for `V`.
"""
TrialFunction(V::FunctionSpace) = Expression(fenics.TrialFunction(V.pyobject))

"""
    TestFunction(V::FunctionSpace)

Construct a symbolic test function for `V`.
"""
TestFunction(V::FunctionSpace) = Expression(fenics.TestFunction(V.pyobject))

"""
    TrialFunctions(V::FunctionSpace)

Return the component trial functions for a mixed function space `V`.
"""
function TrialFunctions(V::FunctionSpace)
    vec = fenics.TrialFunctions(V.pyobject)
    expr_vec = [Expression(elem) for elem in vec]
    return expr_vec
end

"""
    TestFunctions(V::FunctionSpace)

Return the component test functions for a mixed function space `V`.
"""
function TestFunctions(V::FunctionSpace)
    vec = fenics.TestFunctions(V.pyobject)
    expr_vec = [Expression(elem) for elem in vec]
    return expr_vec
end

#Below are attributes for the argument *class*

export Argument, TrialFunction, TrialFunctions, TestFunction, TestFunctions

"""
    Constant

Abstract wrapper for a FEniCS constant expression.
"""
abstract type Constant <: fenicsobject end

struct ConstantImpl <: Constant
    pyobject::PyObject
end

Constant(pyobject::PyObject) = ConstantImpl(pyobject)
"""
    Constant(x::Union{Real, Tuple})

Construct a constant symbolic expression from a scalar or tuple value.

# Arguments

- `x`: Scalar or vector value represented by the constant.
"""
Constant(x::Union{Real, Tuple}) = Expression(fenics.Constant(x, name = "Constant($x)"))
export Constant

"""
    FeFunction

Abstract wrapper for a FEniCS finite-element function.
"""
abstract type FeFunction <: fenicsobject end

struct FeFunctionImpl <: FeFunction
    pyobject::PyObject
end

FeFunction(pyobject::PyObject) = FeFunctionImpl(pyobject)
"""
    FeFunction(V::FunctionSpace; name::String = "")

Construct a finite-element function on `V`.

# Arguments

- `V`: Function space containing the function.

# Keyword Arguments

- `name`: Optional name passed to FEniCS. The default creates an unnamed
  function.
"""
function FeFunction(V::FunctionSpace; name::String = "")
    if name == ""
        return FeFunction(fenics.Function(V.pyobject))
    else
        return FeFunction(fenics.Function(V.pyobject, name = name))
    end
end
"""
    assign(solution1::FeFunction, solution2)
    assign(solution::FeFunction, data::AbstractArray)

Assign values from another finite-element function or an array into a
finite-element function.

# Arguments

- `solution1`, `solution`: Destination finite-element function.
- `solution2`: Source finite-element function.
- `data`: Array of local coefficient values.
"""
function assign(solution1::FeFunction, solution2)
    return fenicspycall(solution1, :assign, solution2.pyobject)
end

function assign(solution::FeFunction, data::AbstractArray)
    return solution.pyobject.vector().set_local(data)
end

"""
    geometric_dimension(expr::Union{FeFunction, Expression})

Return the geometric dimension associated with `expr`.
"""
function geometric_dimension(expr::Union{FeFunction, Expression})
    return fenicspycall(expr, :geometric_dimension)
end
export geometric_dimension

"""
    split(fun::FeFunction)

Split a mixed FEniCS finite-element function into its component functions.

# Arguments
- `fun`: Mixed finite-element function to split.

# Examples
```julia
julia> components = split(mixed_solution);
```
"""
function split(fun::FeFunction)
    vec = fenics.split(fun.pyobject)
    expr_vec = [FeFunction(spl) for spl in vec]
    return expr_vec
end

"""
    py_split(fun::FeFunction)

Split `fun` using the wrapped Python method and return its component
functions.
"""
function py_split(fun::FeFunction)
    vec = fun.pyobject.split()
    expr_vec = [FeFunction(spl) for spl in vec]
    return expr_vec
end

export FeFunction, assign, split, py_split

"""
    Expression(cppcode; kw...)

Construct a symbolic FEniCS expression from C++ expression code.

# Arguments

- `cppcode`: Expression code accepted by FEniCS.
- `kw...`: Keyword arguments forwarded to the FEniCS expression constructor,
  such as `degree`.
"""
function Expression(cppcode; kw...)
    return Expression(fenics.Expression(cppcode; kw...))
end
"""
    Identity(dim::Int)

Construct a symbolic identity tensor of dimension `dim`.
"""
Identity(dim::Int) = Expression(fenics.Identity(dim))
"""
    inner(u, v)

Construct the symbolic inner product of two FEniCS expressions.

# Arguments

- `u`, `v`: Symbolic expressions or finite-element functions.
"""
function inner(u::Union{Expression, FeFunction}, v::Union{Expression, FeFunction})
    return Expression(fenics.inner(u.pyobject, v.pyobject))
end
"""
    outer(u, v)

Construct the symbolic outer product of two FEniCS expressions.
"""
function outer(u::Union{Expression, FeFunction}, v::Union{Expression, FeFunction})
    return Expression(fenics.outer(u.pyobject, v.pyobject))
end
"""
    dot(u, v)

Construct the symbolic dot product of two FEniCS expressions.
"""
function dot(u::Union{Expression, FeFunction}, v::Union{Expression, FeFunction})
    return Expression(fenics.dot(u.pyobject, v.pyobject))
end
"""
    grad(u)
    ∇(u)

Construct the symbolic gradient of `u`.
"""
grad(u::Union{Expression, FeFunction}) = Expression(fenics.grad(u.pyobject))
"""
    ∇(u)

Construct the symbolic gradient of `u`.
"""
∇(u::Union{Expression, FeFunction}) = Expression(fenics.grad(u.pyobject))

"""
    nabla_grad(u)

Construct the UFL nabla gradient of `u`.
"""
nabla_grad(u::Union{Expression, FeFunction}) = Expression(ufl.nabla_grad(u.pyobject))

"""
    nabla_div(u)

Construct the UFL nabla divergence of `u`.
"""
nabla_div(u::Union{Expression, FeFunction}) = Expression(ufl.nabla_div(u.pyobject))
"""
    div(u::Union{Expression, FeFunction})

Construct the FEniCS symbolic divergence of `u`.

# Arguments
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> divergence = div(vector_expression);
```
"""
div(u::Union{Expression, FeFunction}) = Expression(fenics.div(u.pyobject))
"""
    cross(u, v)

Construct the symbolic cross product of two FEniCS expressions.
"""
function cross(u::Union{Expression, FeFunction}, v::Union{Expression, FeFunction})
    return Expression(fenics.cross(u.pyobject, v.pyobject))
end
"""
    tr(u)

Construct the symbolic trace of `u`.
"""
tr(u::Union{Expression, FeFunction}) = Expression(fenics.tr(u.pyobject))
"""
    sqrt(u::Union{Expression, FeFunction})

Construct the symbolic square root of a FEniCS expression or function.

# Arguments
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> magnitude = sqrt(inner(gradient, gradient));
```
"""
sqrt(u::Union{Expression, FeFunction}) = Expression(fenics.sqrt(u.pyobject))
"""
    sym(u)

Construct the symmetric part of `u`.
"""
sym(u::Union{Expression, FeFunction}) = Expression(fenics.sym(u.pyobject))

"""
    len(u)

Return the Python length of a FEniCS expression or finite-element function.
"""
len(U::Union{Expression, FeFunction}) = length(U.pyobject)

sin(u::Union{Expression, FeFunction}) = Expression(fenics.sin(u.pyobject))
cos(u::Union{Expression, FeFunction}) = Expression(fenics.cos(u.pyobject))
tan(u::Union{Expression, FeFunction}) = Expression(fenics.tan(u.pyobject))
asin(u::Union{Expression, FeFunction}) = Expression(fenics.asin(u.pyobject))
acos(u::Union{Expression, FeFunction}) = Expression(fenics.acos(u.pyobject))
atan(u::Union{Expression, FeFunction}) = Expression(fenics.atan(u.pyobject))
exp(u::Union{Expression, FeFunction}) = Expression(fenics.exp(u.pyobject))
log(u::Union{Expression, FeFunction}) = Expression(fenics.ln(u.pyobject))

"""
    besseli(nu::Int, u::Union{Expression, FeFunction})

Construct the modified Bessel function of the first kind for a FEniCS value.

# Arguments
- `nu`: Integer Bessel order.
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> radial_mode = besseli(0, radius);
```
"""
function besseli(nu::Int, u::Union{Expression, FeFunction})
    return Expression(fenics.bessel_I(nu, u.pyobject))
end
"""
    besselj(nu::Int, u::Union{Expression, FeFunction})

Construct the Bessel function of the first kind for a FEniCS value.

# Arguments
- `nu`: Integer Bessel order.
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> radial_mode = besselj(0, radius);
```
"""
function besselj(nu::Int, u::Union{Expression, FeFunction})
    return Expression(fenics.bessel_J(nu, u.pyobject))
end
"""
    besselk(nu::Int, u::Union{Expression, FeFunction})

Construct the modified Bessel function of the second kind for a FEniCS value.

# Arguments
- `nu`: Integer Bessel order.
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> radial_mode = besselk(0, radius);
```
"""
function besselk(nu::Int, u::Union{Expression, FeFunction})
    return Expression(fenics.bessel_K(nu, u.pyobject))
end
"""
    bessely(nu::Int, u::Union{Expression, FeFunction})

Construct the Bessel function of the second kind for a FEniCS value.

# Arguments
- `nu`: Integer Bessel order.
- `u`: FEniCS symbolic expression or finite-element function.

# Examples
```julia
julia> radial_mode = bessely(0, radius);
```
"""
function bessely(nu::Int, u::Union{Expression, FeFunction})
    return Expression(fenics.bessel_Y(nu, u.pyobject))
end

"""
    interpolate(solution1::FeFunction, solution2::Expression)
    interpolate(ex, V::FunctionSpace)

Interpolate an expression into a finite-element function or function space.

# Arguments

- `solution1`: Destination finite-element function for the first method.
- `solution2`: Expression to interpolate for the first method.
- `ex`: Expression to interpolate for the second method.
- `V`: Destination function space for the second method.
"""
function interpolate(solution1::FeFunction, solution2::Expression)
    return FeFunction(fenicspycall(solution1, :interpolate, solution2.pyobject))
end

Expression(x::FEniCS.Expression) = convert(Expression, x)
export Expression, Identity, inner, grad, nabla_grad, nabla_div, div, outer, dot, cross, tr,
    sqrt, sym, len, interpolate
export ∇

"""
    compute_vertex_values(expr, mesh::Mesh)

Compute the values of `expr` at the vertices of `mesh`.

# Arguments

- `expr`: FEniCS expression or finite-element function.
- `mesh`: Mesh whose vertices are used for evaluation.
"""
function compute_vertex_values(expr::Expression, mesh::Mesh)
    return fenicspycall(expr, :compute_vertex_values, mesh.pyobject)
end
function compute_vertex_values(expr::FeFunction, mesh::Mesh)
    return fenicspycall(expr, :compute_vertex_values, mesh.pyobject)
end

export compute_vertex_values

"""
    Measure

Abstract wrapper for a FEniCS integration measure.
"""
abstract type Measure <: fenicsobject end

struct MeasureImpl <: Measure
    pyobject::PyObject
end

Measure(pyobject::PyObject) = MeasureImpl(pyobject)
"""
    directional_derivative(solution1::FeFunction, direction)

Return the directional derivative of `solution1` in `direction`.
"""
function directional_derivative(solution1::FeFunction, direction)
    return FeFunction(fenicspycall(solution1, :dx, direction))
end

"""
    dx

FEniCS cell integration measure used to construct a variational form.

Use `dx` to integrate over the cells of a mesh, for example
`inner(grad(u), grad(v)) * dx`.
"""
dx = nothing

"""
    ds

FEniCS exterior-facet integration measure used to construct a variational
form.

Use `ds` to integrate over the exterior boundary of a mesh.
"""
ds = nothing

"""
    dS

FEniCS interior-facet integration measure used to construct a variational
form.

Use `dS` to integrate over interior facets in a discontinuous Galerkin form.
"""
dS = nothing

"""
    dP

FEniCS point integration measure used to construct a variational form.

Use `dP` when the form integrates over point entities supported by the
underlying FEniCS installation.
"""
dP = nothing

export dx, ds, dS, dP, directional_derivative

#https://github.com/FEniCS/Expression/blob/master/Expression/measure.py
"""
    Form

Abstract wrapper for a FEniCS variational form.
"""
abstract type Form <: fenicsobject end

struct FormImpl <: Form
    pyobject::PyObject
end

Form(pyobject::PyObject) = FormImpl(pyobject)

function *(expr::Union{Expression, FeFunction}, measure::Measure)
    return Expression(measure.pyobject.__rmul__(expr.pyobject))
end

function *(expr::Union{Expression, FeFunction}, expr2::Union{Expression, FeFunction})
    return Expression(expr.pyobject.__mul__(expr2.pyobject))
end
function *(expr::Real, expr2::Union{Expression, FeFunction})
    return Expression(expr2.pyobject.__mul__(expr))
end
function *(expr::Union{Expression, FeFunction}, expr2::Real)
    return Expression(expr.pyobject.__mul__(expr2))
end

function +(expr::Union{Expression, FeFunction}, expr2::Real)
    return Expression(expr.pyobject.__add__(expr2))
end
function +(expr::Real, expr2::Union{Expression, FeFunction})
    return Expression(expr2.pyobject.__add__(expr))
end
function +(expr::Union{Expression, FeFunction}, expr2::Union{Expression, FeFunction})
    return Expression(expr.pyobject.__add__(expr2.pyobject))
end

function -(expr::Union{Expression, FeFunction}, expr2::Real)
    return Expression(expr.pyobject.__sub__(expr2))
end
function -(expr::Real, expr2::Union{Expression, FeFunction})
    return -1 * (Expression(expr2.pyobject.__sub__(expr)))
end
function -(expr::Union{Expression, FeFunction}, expr2::Union{Expression, FeFunction})
    return Expression(expr.pyobject.__sub__(expr2.pyobject))
end

function /(expr::Union{Expression, FeFunction}, expr2::Real)
    return Expression(expr.pyobject.__div__(expr2))
end
function /(expr::Union{Expression, FeFunction}, expr2::Union{Expression, FeFunction})
    return Expression(expr.pyobject.__div__(expr2.pyobject))
end

function /(expr::Real, expr2::Union{Expression, FeFunction})
    x = expr2 * expr2
    y = x / expr
    z = expr2 / y
    return Expression(z)
end

function ^(expr::Union{Expression, FeFunction}, expr2::Real)
    return Expression(expr.pyobject.__pow__(expr2))
end
function ^(expr::Union{Expression, FeFunction}, expr2::Union{Expression, FeFunction})
    return Expression(expr.pyobject.__pow__(expr2.pyobject))
end

"""
    Transpose(object::Expression)

Return the transpose of a symbolic FEniCS expression.

# Arguments

- `object`: Expression to transpose.
"""
function Transpose(object::Expression)
    x = object.pyobject.T
    y = Expression(x)
    return y
end
export Transpose

"""
    rhs(equation::Expression)

Extract the right-hand side from a combined bilinear and linear form. The
linear part is negated by the FEniCS convention.

# Example

```julia
a = u * v * dx + f * v * dx
L = rhs(a)
```
"""
rhs(equation::Expression) = Expression(fenics.rhs(equation.pyobject))
"""
    lhs(equation::Expression)

Extract the bilinear left-hand side from a combined variational form.

# Example

```julia
a = u * v * dx + f * v * dx
A = lhs(a)
```
"""
lhs(equation::Expression) = Expression(fenics.lhs(equation.pyobject))

export lhs, rhs
#this assembles the matrix from a fenics form
"""
    Matrix

Abstract wrapper for a FEniCS assembled matrix.
"""
abstract type Matrix <: fenicsobject end

struct MatrixImpl <: Matrix
    pyobject::PyObject
end

Matrix(pyobject::PyObject) = MatrixImpl(pyobject)

"""
    Matrix(a::T) where {T <: Real}

Return a scalar unchanged when a FEniCS form reduces to a real number.
"""
Matrix(a::T) where {T <: Real} = a

"""
    assemble(assembly_item; tensor = nothing,
        form_compiler_parameters = nothing, add_values = false,
        finalize_tensor = true, keep_diagonal = false, backend = nothing)

Assemble a FEniCS form or expression into a matrix-like object.

# Arguments

- `assembly_item`: Form or expression to assemble.

# Keyword Arguments

- `tensor`: Optional existing tensor to fill.
- `form_compiler_parameters`: Optional FEniCS compiler parameters.
- `add_values`: Whether to add into an existing tensor.
- `finalize_tensor`: Whether to finalize the assembled tensor.
- `keep_diagonal`: Whether to preserve the tensor diagonal.
- `backend`: Optional assembly backend.
"""
function assemble(
        assembly_item::Union{Form, Expression}; tensor = nothing,
        form_compiler_parameters = nothing, add_values = false,
        finalize_tensor = true, keep_diagonal = false, backend = nothing
    )
    return Matrix(
        fenics.assemble(
            assembly_item.pyobject,
            tensor = tensor,
            form_compiler_parameters = form_compiler_parameters,
            add_values = add_values, finalize_tensor = finalize_tensor,
            keep_diagonal = keep_diagonal, backend = backend
        )
    )
end
export assemble

"""
    assemble_local(assembly_item::Union{Form, Expression}, cell::Cell)

Assemble `assembly_item` locally on `cell`.

# Arguments

- `assembly_item`: Form or expression to assemble.
- `cell`: Cell on which local assembly is performed.
"""
function assemble_local(assembly_item::Union{Form, Expression}, cell::Cell)
    return fenics.assemble_local(assembly_item.pyobject, cell.pyobject)
end
export assemble_local

#I have changed this to Function+Form

#https://fenicsproject.org/olddocs/dolfin/1.6.0/python/programmers-reference/cpp/fem/DirichletBC.html
"""
    sub_domain

Abstract wrapper for a FEniCS boundary subdomain.
"""
abstract type sub_domain <: fenicsobject end

struct sub_domainImpl <: sub_domain
    pyobject::PyObject
end

sub_domain(pyobject::PyObject) = sub_domainImpl(pyobject)

"""
    BoundaryCondition

Abstract wrapper for a FEniCS boundary condition.
"""
abstract type BoundaryCondition <: fenicsobject end

struct BoundaryConditionImpl <: BoundaryCondition
    pyobject::PyObject
end

BoundaryCondition(pyobject::PyObject) = BoundaryConditionImpl(pyobject)
"""
    DirichletBC(V::FunctionSpace, g, sub_domain)

Construct a Dirichlet boundary condition on `V`.

The boundary value `g` may be a FEniCS expression, number, or tuple. The
`sub_domain` argument identifies the boundary on which the condition applies.

# Arguments

- `V`: Function space constrained by the boundary condition.
- `g`: Boundary value.
- `sub_domain`: Boundary selector accepted by FEniCS.
"""
function DirichletBC(V::FunctionSpace, g, sub_domain)
    return BoundaryCondition(fenics.DirichletBC(V.pyobject, g.pyobject, sub_domain))
end #look this up with example also removed type from g(Should be expression)
function DirichletBC(V::FunctionSpace, g::Number, sub_domain)
    return BoundaryCondition(fenics.DirichletBC(V.pyobject, g, sub_domain))
end #look this up with example also removed type from g(Should be expression)
function DirichletBC(V::FunctionSpace, g::Tuple, sub_domain)
    return BoundaryCondition(fenics.DirichletBC(V.pyobject, g, sub_domain))
end #look this up with example also removed type from g(Should be expression)

export DirichletBC
"""
    CompiledSubDomain(cppcode::String)

Compile a C++ boundary predicate into a FEniCS subdomain object.

# Arguments

- `cppcode`: C++ predicate source accepted by FEniCS.
"""
CompiledSubDomain(cppcode::String) = sub_domain(fenics.CompiledSubDomain(cppcode))
export CompiledSubDomain

"""
    apply(bcs::BoundaryCondition, matrix::Matrix)

Apply a boundary condition to an assembled matrix.

# Arguments

- `bcs`: Boundary condition to apply.
- `matrix`: Assembled matrix to modify.
"""
apply(bcs::BoundaryCondition, matrix::Matrix) = fenicspycall(bcs, :apply, matrix.pyobject)
apply(bcs, matrix::Matrix) = fenicspycall(BoundaryCondition(bcs), :apply, matrix.pyobject)

export apply

"""
    assemble_system_julia(a::Expression, L::Expression, bc = nothing)

Assemble a variational system and return Julia arrays.

# Arguments

- `a`: Bilinear form.
- `L`: Linear form.
- `bc`: Optional boundary condition or collection of boundary conditions.

# Returns

A tuple `(A, b)` containing the assembled Julia matrix and right-hand side.
"""
function assemble_system_julia(a::Expression, L::Expression)
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject)
    A = A_fenics[:array]()
    b = b_fenics[:array]()
    return A, b
end

function assemble_system_julia(a::Expression, L::Expression, bc)
    bcs_py = [bcs.pyobject for bcs in bc]
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject, bcs_py)
    A = A_fenics[:array]()
    b = b_fenics[:array]()
    return A, b
end

function assemble_system_julia(a::Expression, L::Expression, bc::BoundaryCondition)
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject, bc.pyobject)
    A = A_fenics[:array]()
    b = b_fenics[:array]()
    return A, b
end

##add assemble_system to FEniCS objects for solving

"""
    assemble_system(a::Expression, L::Expression, bc = nothing)

Assemble a variational system and return FEniCS-backed matrix objects.

# Arguments

- `a`: Bilinear form.
- `L`: Linear form.
- `bc`: Optional boundary condition or collection of boundary conditions.

# Returns

A tuple `(A, b)` containing the assembled system matrix and right-hand side.
"""
function assemble_system(a::Expression, L::Expression)
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject)
    A = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A, b
end

function assemble_system(a::Expression, L::Expression, bc)
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject, bc)
    A = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A, b
end

function assemble_system(a::Expression, L::Expression, bc::BoundaryCondition)
    A_fenics, b_fenics = fenics.assemble_system(a.pyobject, L.pyobject, bc)
    A = Matrix(A_fenics)
    b = Matrix(b_fenics)
    return A, b
end

export assemble_system, assemble_system_julia

"""
 For a full list of supported arguments, and their usage
please refer to http://matplotlib.org/api/pyplot_api.html
not all kwargs have been imported. Should you require any that are not imported
open as issue, and I will attempt to add them.
Deprecate this in a future version
"""
function Plot(
        in_plot::Union{Mesh, FunctionSpace, FeFunction}; alpha = 1, animated = false,
        antialiased = true, color = "grey", dash_capstyle = "butt",
        dash_joinstyle = "miter", dashes = "", drawstyle = "default",
        fillstyle = "full", label = "s", linestyle = "solid", linewidth = 1,
        marker = "", markeredgecolor = "grey", markeredgewidth = "",
        markerfacecolor = "grey", markerfacecoloralt = "grey", markersize = 1,
        markevery = "none", visible = true, title = ""
    )
    return fenics.common.plotting.plot(
        in_plot.pyobject,
        alpha = alpha, animated = animated,
        antialiased = antialiased, color = color,
        dash_capstyle = dash_capstyle,
        dash_joinstyle = dash_joinstyle, dashes = dashes,
        drawstyle = drawstyle, fillstyle = fillstyle, label = label,
        linestyle = linestyle, linewidth = linewidth,
        marker = marker, markeredgecolor = markeredgecolor,
        markeredgewidth = markeredgewidth,
        markerfacecolor = markerfacecolor,
        markerfacecoloralt = markerfacecoloralt,
        markersize = markersize, markevery = markevery,
        visible = visible, title = title
    ) #the first is the keyword argument, the second is the value
end #the first is the keyword argument, the second is the value
#export Plot

"""
    FiniteElement

Abstract wrapper for a FEniCS finite element.
"""
abstract type FiniteElement <: fenicsobject end

struct FiniteElementImpl <: FiniteElement
    pyobject::PyObject
end

FiniteElement(pyobject::PyObject) = FiniteElementImpl(pyobject)
"""
    FiniteElement(family::StringOrSymbol, cell = nothing, degree = nothing,
        form_degree = nothing, quad_scheme = nothing, variant = nothing)

Construct a FEniCS finite element.

# Arguments

- `family`: Finite-element family name.
- `cell`: Geometric cell, such as [`triangle`](@ref).
- `degree`: Polynomial degree.
- `form_degree`: Optional FEEC form degree.
- `quad_scheme`: Optional quadrature scheme.
- `variant`: Optional local-basis variant.
"""
function FiniteElement(
        family::StringOrSymbol, cell = nothing, degree = nothing,
        form_degree = nothing, quad_scheme = nothing, variant = nothing
    )
    return FiniteElement(
        fenics.FiniteElement(
            family = family, cell = cell, degree = degree,
            form_degree = form_degree, quad_scheme = quad_scheme,
            variant = variant
        )
    )
end
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

"""
    hexahedron

FEniCS cell object describing a hexahedral finite-element cell.

Use it as the `cell` argument to [`FiniteElement`](@ref) when constructing a
hexahedral element.
"""
hexahedron = nothing

"""
    tetrahedron

FEniCS cell object describing a tetrahedral finite-element cell.

Use it as the `cell` argument to [`FiniteElement`](@ref) when constructing a
tetrahedral element.
"""
tetrahedron = nothing

"""
    quadrilateral

FEniCS cell object describing a quadrilateral finite-element cell.

Use it as the `cell` argument to [`FiniteElement`](@ref) when constructing a
quadrilateral element.
"""
quadrilateral = nothing

"""
    triangle

FEniCS cell object describing a triangular finite-element cell.

Use it as the `cell` argument to [`FiniteElement`](@ref) when constructing a
triangular element.
"""
triangle = nothing

export hexahedron, tetrahedron, quadrilateral, triangle

"""
    family(finiteelement::FiniteElement)

Return the family name of `finiteelement`.
"""
family(finiteelement::FiniteElement) = fenicspycall(finiteelement, :family)

"""
    cell(finiteelement::FiniteElement)

Return the geometric cell of `finiteelement`.
"""
cell(finiteelement::FiniteElement) = fenicspycall(finiteelement, :cell)

"""
    degree(finiteelement::FiniteElement)

Return the polynomial degree of `finiteelement`.
"""
degree(finiteelement::FiniteElement) = fenicspycall(finiteelement, :degree)

#form_degree(finiteelement::FiniteElement) = fenicspycall(finiteelement, :form_degree)

#quad_scheme(finiteelement::FiniteElement) = fenicspycall(finiteelement, :quad_scheme)

"""
    variant(finiteelement::FiniteElement)

Return the local-basis variant of `finiteelement`.
"""
variant(finiteelement::FiniteElement) = fenicspycall(finiteelement, :variant)

"""
    reconstruct(finiteelement::FiniteElement; family = nothing,
        cell = nothing, degree = nothing)

Construct a new finite element with selected properties replaced.

# Keyword Arguments

- `family`: Replacement family, or `nothing` to preserve the current family.
- `cell`: Replacement cell, or `nothing` to preserve the current cell.
- `degree`: Replacement degree, or `nothing` to preserve the current degree.
"""
function reconstruct(
        finiteelement::FiniteElement, ; family = nothing, cell = nothing,
        degree = nothing
    )
    return FiniteElement(fenicspycall(finiteelement, :reconstruct, family, cell, degree))
end

"""
    sobolev_space(finiteelement::FiniteElement)

Return the Sobolev space associated with `finiteelement`.
"""
sobolev_space(finiteelement::FiniteElement) = fenicspycall(finiteelement, :sobolev_space)
export family, cell, degree, variant, reconstruct, sobolev_space

"""
    FacetNormal(mesh::Mesh)

Construct the outward facet-normal expression for `mesh`.
"""
FacetNormal(mesh::Mesh) = Expression(fenics.FacetNormal(mesh.pyobject))

export FacetNormal

"""
    MixedElement

Abstract wrapper for a FEniCS mixed finite element.
"""
abstract type MixedElement <: fenicsobject end

struct MixedElementImpl <: MixedElement
    pyobject::PyObject
end

MixedElement(pyobject::PyObject) = MixedElementImpl(pyobject)

"""
    MixedElement(vec::Array{FEniCS.FiniteElementImpl, 1})

Construct a mixed finite element from a vector of finite-element wrappers.

# Arguments

- `vec`: Finite-element components to combine.
"""
function MixedElement(vec::Array{FEniCS.FiniteElementImpl, 1})
    pyvec = [elem.pyobject for elem in vec]
    return me = MixedElement(fenics.MixedElement(pyvec))
end

export MixedElement

"""
    FunctionSpace(mesh::Mesh, element::Union{FiniteElement, MixedElement})

Construct a function space from a finite element or mixed finite element.
"""
function FunctionSpace(mesh::Mesh, element::Union{FiniteElement, MixedElement})
    return FunctionSpace(fenics.FunctionSpace(mesh.pyobject, element.pyobject))
end

export FunctionSpace
