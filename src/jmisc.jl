#this file contains miscallaneous functions mainly related to the solve.jl file.

#   https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/fem/solving/solve.html

"""
lists available lu solver methods
"""
list_lu_solver_methods()=fenics.list_lu_solver_methods()
"""
lists available krylov solver methods
"""
list_krylov_solver_methods()=fenics.list_krylov_solver_methods()
"""
lists available krylov solver preconditioners
"""
list_krylov_solver_preconditioners()=fenics.list_krylov_solver_preconditioners()
"""
lists available linear solver methods
"""
list_linear_solver_methods()=fenics.list_linear_solver_methods()


"""
Provide values for some of the constants
"""
#not sure if this was a necessary import, but maybe consistency is required.
DOLFIN_PI()=fenics.DOLFIN_PI
DOLFIN_EPS()=fenics.DOLFIN_EPS
DOLFIN_SQRT_EPS()=fenics.DOLFIN_SQRT_EPS
