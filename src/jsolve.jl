#this file contains functions/wrappers related to the solve function in FEniCS
#https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/fem/solving/solve.html
using FEniCS

solve(A::Matrix,x, b::Matrix,solvers...) =fenics.solve(A.pyobject,x,b.pyobject,solvers...)
export solve

#lvsolve is the linear variational solver
lvsolve(a,L,u;solver_parameters::Dict=Dict("linear_solver"=>"default"),form_compiler_parameters::Dict=Dict("optimize"=>true))=fenics.solve(a.pyobject==L.pyobject, u.pyobject,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)

lvsolve(a,L,u,bcs=nothing;solver_parameters::Dict=Dict("linear_solver"=>"default"),form_compiler_parameters::Dict=Dict("optimize"=>true))=fenics.solve(a.pyobject==L.pyobject, u.pyobject, bcs=bcs.pyobject,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)
#allows BoundaryCondition to be provided in an AbstractArray (of type BoundaryCondition)
function lvsolve(a,L,u,bcs::AbstractArray;solver_parameters::Dict=Dict("linear_solver"=>"default"),form_compiler_parameters::Dict=Dict("optimize"=>true))
    bcs_py = [bc.pyobject for bc in bcs]
    fenics.solve(a.pyobject==L.pyobject, u.pyobject, bcs=bcs_py,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)
end

export lvsolve
#Dict("linear_solver"=>"default")
#Dict("optimize"=>true)
#nlvsolve is the non-linear variational solver
nlvsolve(F,u;J=nothing,solver_parameters::Dict=Dict("nonlinear_solver"=>"newton"),form_compiler_parameters::Dict=Dict("optimize"=>true))=fenics.solve(F.pyobject==0,u.pyobject,J=J,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)
nlvsolve(F,u,bcs=nothing;J=nothing,solver_parameters::Dict=Dict("nonlinear_solver"=>"newton"),form_compiler_parameters::Dict=Dict("optimize"=>true))=fenics.solve(F.pyobject==0,u.pyobject,J=J,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)

nlvsolve(F,u,bcs::BoundaryCondition;J=nothing,solver_parameters::Dict=Dict("nonlinear_solver"=>"newton"),form_compiler_parameters::Dict=Dict("optimize"=>true))=fenics.solve(F.pyobject==0,u.pyobject,bcs=bcs.pyobject,J=J,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)
#allows BoundaryCondition to be provided in an AbstractArray (of type BoundaryCondition)
function nlvsolve(F,u,bcs::AbstractArray;J=nothing,solver_parameters::Dict=Dict("nonlinear_solver"=>"newton"),form_compiler_parameters::Dict=Dict("optimize"=>true))
    bcs_py = [bc.pyobject for bc in bcs]
    fenics.solve(F.pyobject==0,u.pyobject,bcs=bcs,J=J,solver_parameters=solver_parameters,form_compiler_parameters=form_compiler_parameters)
end



export nlvsolve

#anlvsolve corresponds to the adapative non-linear solver.
anlvsolve(F,a,u,bcs,tol,M)=fenics.solve(F.pyobject==a.pyobject,u.pyobject,bcs=bcs.pyobject,tol=tol,M=M)
#this function hasnt been tested yet, so isnt exported



"""
errornorm is the function to calculate the error between our exact and calculated
solution. The norm kwarg defines the norm measure used (by default it is the L2 norm)
"""
errornorm(ans,sol;norm="L2") = fenics.errornorm(ans.pyobject,sol.pyobject,norm)
export errornorm

File(path::StringOrSymbol)=fenics.File(path) #used to store the solution in various formats

function File(path::StringOrSymbol,object::FeFunction)
    vtkfile = File(path)
    vtkfile << object.pyobject
end

function File(path::StringOrSymbol,object::FeFunction, time::Number)
    vtkfile = File(path)
    vtkfile << (object.pyobject,time)
end


function File(path::StringOrSymbol,object::Mesh)
    vtkfile = File(path)
    vtkfile << object.pyobject
end

function File(path::StringOrSymbol,object::Mesh, time::Number)
    vtkfile = File(path)
    vtkfile << (object.pyobject,time)
end


export File



XDMFFile(path::StringOrSymbol) = fenics.XDMFFile(path)
export XDMFFile

TimeSeries(path::StringOrSymbol) = fenics.TimeSeries(path)
retrieve(timeseries,placeholder,time) = timeseries.retrieve(placeholder,time)
export TimeSeries, retrieve


write(path,solution,time::Number) = path.write(solution.pyobject, time)
write(path,solution::PyObject,time::Number) = path.write(solution, time)

store(path,solution,time::Number) = path.store(solution.pyobject, time)
store(path,solution::PyObject,time::Number) = path.store(solution, time)

export write,store

array(matrix) = fenicspycall(matrix, :gather_on_zero)
vector(solution::FeFunction) = fenicspycall(solution,:vector) #
interpolate(ex, V::FunctionSpace) = FeFunction(fenics.interpolate(ex.pyobject, V.pyobject))

export vector, interpolate,array


"the following function is used to extract the array from a solution/form"
function get_array(form::Expression)
    assembled_form = assemble(form)
    return array(assembled_form)
end

function get_array(solution::FeFunction)
    generic_vector = vector(solution)
    instantiated_vector = fenics.Vector(generic_vector)
    return instantiated_vector.gather_on_zero()
end
"""
we return the array from an assembled form
"""
function get_array(assembled_form::Matrix)
    return array(assembled_form)
end

export get_array

"""
Return projection of given expression *v* onto the finite element space *V*
*Example of usage*
          v = Expression("sin(pi*x[0])")
          V = FunctionSpace(mesh, "Lagrange", 1)
          Pv = project(v, V)
"""
project(v::Union{FeFunction,Expression},V::FunctionSpace)  = FeFunction(fenics.project(v.pyobject,V.pyobject))
export project
