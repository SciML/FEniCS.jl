#this file contains functions/wrappers related to the solve function in FEniCS
#https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/fem/solving/solve.html
using FEniCS

solve(A::Matrix,x, b::Matrix,solvers...) =fenics.solve(A.pyobject,x,b.pyobject,solvers...)
export solve

#lvsolve is the linear variational solver
lvsolve(a,L,u)=fenics.solve(a.pyobject==L.pyobject, u.pyobject)
lvsolve(a,L,u,bcs)=fenics.solve(a.pyobject==L.pyobject, u.pyobject, bcs=bcs.pyobject)
export lvsolve

#nlvsolve is the non-linear variational solver
nlvsolve(F,u;bcs=nothing)=fenics.solve(F.pyobject==0,u.pyobject,bcs=bcs.pyobject)
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

File(path::String)=fenics.File(path) #used to store the solution in various formats
export File

XDMFFile(path::String) = fenics.XDMFFile(path)
export XDMFFile

TimeSeries(path::String) = fenics.TimeSeries(path)
export TimeSeries

array(matrix) = fenicspycall(matrix, :array)
vector(solution) = fenicspycall(solution,:vector) #
Vector(solution) = fenics.Vector(solution) # genericvector fenics
#various overloading of interpolate
interpolate(ex, V::FunctionSpace) = Function(fenics.interpolate(ex.pyobject, V.pyobject))
#interpolate(fun::Function, expr::Expression) = Function(fenics.interpolate(fun.pyobject, expr.pyobject))

export Vector, vector, interpolate,array


"the following function is used to extract the array from a solution/form"
function get_array(form::Form)
    assembled_form = assemble(form)
    return array(assembled_form)
end

function get_array(solution::Function)
    generic_vector = vector(solution)
    instantiated_vector = Vector(generic_vector)
    return instantiated_vector[:array]()
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
project(v::Union{Function,Expression},V::FunctionSpace)  = Function(fenics.project(v.pyobject,V.pyobject))
export project
