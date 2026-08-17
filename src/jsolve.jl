#this file contains functions/wrappers related to the solve function in FEniCS
#https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/fem/solving/solve.html

"""
    solve(A::Matrix, x, b::Matrix, solvers...)

Solve the linear algebraic system represented by `A`, `x`, and `b` using
the wrapped FEniCS solver.

# Arguments

- `A`: Assembled system matrix.
- `x`: Solution vector or FEniCS solution object.
- `b`: Right-hand-side matrix or vector.
- `solvers...`: Additional solver arguments forwarded to FEniCS.
"""
function solve(A::Matrix, x, b::Matrix, solvers...)
    return fenics.solve(A.pyobject, x, b.pyobject, solvers...)
end
export solve

"""
    lvsolve(a, L, u, bcs = nothing;
        solver_parameters = Dict("linear_solver" => "default"),
        form_compiler_parameters = Dict("optimize" => true))

Solve a linear variational problem.

# Arguments

- `a`: Bilinear form.
- `L`: Linear form.
- `u`: Unknown finite-element function.
- `bcs`: Optional boundary condition or collection of boundary conditions.

# Keyword Arguments

- `solver_parameters`: FEniCS linear-solver parameters.
- `form_compiler_parameters`: FEniCS form-compiler parameters.
"""
function lvsolve(
        a, L, u, bcs = nothing;
        solver_parameters::Dict = Dict("linear_solver" => "default"),
        form_compiler_parameters::Dict = Dict("optimize" => true)
    )
    return if bcs === nothing
        fenics.solve(
            a.pyobject == L.pyobject, u.pyobject,
            solver_parameters = solver_parameters,
            form_compiler_parameters = form_compiler_parameters
        )
    else
        fenics.solve(
            a.pyobject == L.pyobject, u.pyobject, bcs = bcs.pyobject,
            solver_parameters = solver_parameters,
            form_compiler_parameters = form_compiler_parameters
        )
    end
end
#allows BoundaryCondition to be provided in an AbstractArray (of type BoundaryCondition)
function lvsolve(
        a, L, u, bcs::AbstractArray;
        solver_parameters::Dict = Dict("linear_solver" => "default"),
        form_compiler_parameters::Dict = Dict("optimize" => true)
    )
    bcs_py = [bc.pyobject for bc in bcs]
    return fenics.solve(
        a.pyobject == L.pyobject, u.pyobject, bcs = bcs_py,
        solver_parameters = solver_parameters,
        form_compiler_parameters = form_compiler_parameters
    )
end

export lvsolve
#Dict("linear_solver"=>"default")
#Dict("optimize"=>true)
"""
    nlvsolve(F, u, bcs = nothing; J = nothing,
        solver_parameters = Dict("nonlinear_solver" => "newton"),
        form_compiler_parameters = Dict("optimize" => true))

Solve a nonlinear variational problem.

# Arguments

- `F`: Nonlinear residual form.
- `u`: Unknown finite-element function.
- `bcs`: Optional boundary condition or collection of boundary conditions.
- `J`: Optional Jacobian form.

# Keyword Arguments

- `solver_parameters`: FEniCS nonlinear-solver parameters.
- `form_compiler_parameters`: FEniCS form-compiler parameters.
"""
function nlvsolve(
        F, u, bcs = nothing; J = nothing,
        solver_parameters::Dict = Dict("nonlinear_solver" => "newton"),
        form_compiler_parameters::Dict = Dict("optimize" => true)
    )
    return if bcs === nothing
        fenics.solve(
            F.pyobject == 0, u.pyobject, J = J, solver_parameters = solver_parameters,
            form_compiler_parameters = form_compiler_parameters
        )
    else
        fenics.solve(
            F.pyobject == 0, u.pyobject, bcs = bcs.pyobject, J = J,
            solver_parameters = solver_parameters,
            form_compiler_parameters = form_compiler_parameters
        )
    end
end

function nlvsolve(
        F, u, bcs::BoundaryCondition; J = nothing,
        solver_parameters::Dict = Dict("nonlinear_solver" => "newton"),
        form_compiler_parameters::Dict = Dict("optimize" => true)
    )
    return fenics.solve(
        F.pyobject == 0, u.pyobject, bcs = bcs.pyobject, J = J,
        solver_parameters = solver_parameters,
        form_compiler_parameters = form_compiler_parameters
    )
end
#allows BoundaryCondition to be provided in an AbstractArray (of type BoundaryCondition)
function nlvsolve(
        F, u, bcs::AbstractArray; J = nothing,
        solver_parameters::Dict = Dict("nonlinear_solver" => "newton"),
        form_compiler_parameters::Dict = Dict("optimize" => true)
    )
    bcs_py = [bc.pyobject for bc in bcs]
    return fenics.solve(
        F.pyobject == 0, u.pyobject, bcs = bcs, J = J,
        solver_parameters = solver_parameters,
        form_compiler_parameters = form_compiler_parameters
    )
end

export nlvsolve

"""
    anlvsolve(F, a, u, bcs, tol, M)

Solve an adaptive nonlinear variational problem through the wrapped FEniCS
solver.

This lower-level helper is kept for compatibility and is not exported.
"""
function anlvsolve(F, a, u, bcs, tol, M)
    return fenics.solve(F.pyobject == a.pyobject, u.pyobject, bcs = bcs.pyobject, tol = tol, M = M)
end
#this function hasnt been tested yet, so isnt exported

"""
    norm(u::FeFunction; normType = "L2", mesh = nothing)

Compute a FEniCS norm of a finite-element function.

# Arguments
- `u`: Finite-element function whose norm is computed.

# Keyword Arguments
- `normType = "L2"`: FEniCS norm identifier.
- `mesh = nothing`: Optional mesh used by FEniCS for the norm computation.

# Examples
```julia
julia> norm(u; normType = "H1")
1.0
```
"""
function norm(u::FeFunction; normType = "L2", mesh::Union{Nothing, Mesh} = nothing)
    if isa(mesh, Nothing)
        return fenics.norm(u.pyobject, normType)
    else
        return fenics.norm(u.pyobject, normType, mesh.pyobject)
    end
end

"""
    errornorm(ans, sol; norm = "L2")

Compute the FEniCS error norm between an exact solution `ans` and a computed
solution `sol`.

# Keyword Arguments

- `norm`: FEniCS norm identifier. The default is `"L2"`.
"""
errornorm(ans, sol; norm = "L2") = fenics.errornorm(ans.pyobject, sol.pyobject, norm)
export errornorm

"""
    File(path::StringOrSymbol)

Create a FEniCS output file at `path`.
"""
File(path::StringOrSymbol) = fenics.File(path) #used to store the solution in various formats

function File(path::StringOrSymbol, object::FeFunction)
    vtkfile = File(path)
    return vtkfile << object.pyobject
end

function File(path::StringOrSymbol, object::FeFunction, time::Number)
    vtkfile = File(path)
    return vtkfile << (object.pyobject, time)
end

function File(path::StringOrSymbol, object::Mesh)
    vtkfile = File(path)
    return vtkfile << object.pyobject
end

function File(path::StringOrSymbol, object::Mesh, time::Number)
    vtkfile = File(path)
    return vtkfile << (object.pyobject, time)
end

export File

"""
    XDMFFile(path::StringOrSymbol)

Create a FEniCS XDMF output object at `path`.
"""
XDMFFile(path::StringOrSymbol) = fenics.XDMFFile(path)
export XDMFFile

"""
    TimeSeries(path::StringOrSymbol)

Create a FEniCS time-series storage object at `path`.
"""
TimeSeries(path::StringOrSymbol) = fenics.TimeSeries(path)

"""
    retrieve(timeseries, placeholder, time)

Retrieve the value associated with `placeholder` at `time` from a FEniCS
time series.
"""
retrieve(timeseries, placeholder, time) = timeseries.retrieve(placeholder, time)
export TimeSeries, retrieve

"""
    write(path::PyObject, solution::fenicsobject, time::Number)

Write a FEniCS mesh or function to an FEniCS `XDMFFile` or `TimeSeries` at `time`.

# Arguments
- `path`: Python-backed FEniCS output object returned by `XDMFFile` or `TimeSeries`.
- `solution`: FEniCS mesh or function to write.
- `time`: Time associated with the output sample.

# Examples
```julia
julia> write(XDMFFile("solution.xdmf"), u, 0.0)
```
"""
write(path::PyObject, solution::fenicsobject, time::Number) = path.write(solution.pyobject, time)

"""
    store(path::PyObject, solution, time::Number)

Store `solution` at `time` in a FEniCS time-series object.
"""
store(path::PyObject, solution, time::Number) = path.store(solution.pyobject, time)
store(path::PyObject, solution::PyObject, time::Number) = path.store(solution, time)

export write, store

"""
    array(matrix)

Gather a FEniCS matrix or vector-like object into a Julia array on rank zero.
"""
array(matrix) = fenicspycall(matrix, :gather_on_zero)

"""
    vector(solution::FeFunction)

Return the vector backing `solution`.
"""
vector(solution::FeFunction) = fenicspycall(solution, :vector) #

export vector, interpolate, array

"""
    get_array(form::Expression)
    get_array(solution::FeFunction)
    get_array(assembled_form::Matrix)

Extract a Julia array from a FEniCS expression, finite-element function, or
assembled matrix.
"""
function get_array(form::Expression)
    assembled_form = assemble(form)
    return array(assembled_form)
end

function get_array(solution::FeFunction)
    generic_vector = vector(solution)
    instantiated_vector = fenics.Vector(generic_vector)
    return instantiated_vector.gather_on_zero()
end
function get_array(assembled_form::Matrix)
    return array(assembled_form)
end

export get_array

"""
    project(v::Union{FeFunction, Expression}, V::FunctionSpace)

Project `v` onto the finite-element space `V`.

# Example

```julia
v = Expression("sin(pi*x[0])", degree = 2)
V = FunctionSpace(mesh, "Lagrange", 1)
Pv = project(v, V)
```
"""
function project(v::Union{FeFunction, Expression}, V::FunctionSpace)
    return FeFunction(fenics.project(v.pyobject, V.pyobject))
end
export project
