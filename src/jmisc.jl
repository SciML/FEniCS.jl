#this file contains miscallaneous functions mainly related to the solve.jl file.

#   https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/fem/solving/solve.html

"""
    list_lu_solver_methods()

Return the LU solver methods available in the FEniCS installation.
"""
list_lu_solver_methods() = fenics.list_lu_solver_methods()
"""
    list_krylov_solver_methods()

Return the Krylov solver methods available in the FEniCS installation.
"""
list_krylov_solver_methods() = fenics.list_krylov_solver_methods()
"""
    list_krylov_solver_preconditioners()

Return the Krylov preconditioners available in the FEniCS installation.
"""
list_krylov_solver_preconditioners() = fenics.list_krylov_solver_preconditioners()
"""
    list_linear_solver_methods()

Return the linear solver methods available in the FEniCS installation.
"""
list_linear_solver_methods() = fenics.list_linear_solver_methods()

"""
    info_NonLinearVariationalSolver()

Print the default parameters for the FEniCS nonlinear variational solver.
"""
function info_NonLinearVariationalSolver()
    return fenics.info(fenics.NonlinearVariationalSolver.default_parameters(), true)
end
"""
    info_LinearVariatonalSolver()

Print the default parameters for the FEniCS linear variational solver.
"""
function info_LinearVariatonalSolver()
    return fenics.info(fenics.LinearVariationalSolver.default_parameters(), true)
end

export list_lu_solver_methods, list_krylov_solver_methods,
    list_krylov_solver_preconditioners, list_linear_solver_methods,
    info_NonLinearVariationalSolver, info_LinearVariatonalSolver

"""
    DOLFIN_PI()

Return the DOLFIN value of pi.
"""
DOLFIN_PI() = fenics.DOLFIN_PI

"""
    DOLFIN_EPS()

Return the DOLFIN floating-point epsilon.
"""
DOLFIN_EPS() = fenics.DOLFIN_EPS

"""
    DOLFIN_SQRT_EPS()

Return the square root of the DOLFIN floating-point epsilon.
"""
DOLFIN_SQRT_EPS() = fenics.DOLFIN_SQRT_EPS

export DOLFIN_PI, DOLFIN_EPS, DOLFIN_SQRT_EPS

"""
    MPI_Comm

Wrapper for a FEniCS MPI communicator object.
"""
struct MPI_Comm <: fenicsobject
    pyobject::PyObject
end

"""
    mpi_comm_world()

Return the FEniCS world MPI communicator.
"""
mpi_comm_world() = MPI_Comm(fenics.mpi_comm_world())

"""
    mpi_comm_self()

Return the FEniCS self MPI communicator.
"""
mpi_comm_self() = MPI_Comm(fenics.mpi_comm_world())

export MPI_Comm, mpi_comm_world, mpi_comm_self
