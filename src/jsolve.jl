#this file contains functions/wrappers related to the solve function in FEniCS
#https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/fem/solving/solve.html
using FEniCS
lvsolve(a,L,u)=fenics.solve(a.pyobject==L.pyobject, u.pyobject)
lvsolve(a,L,u,bcs)=fenics.solve(a.pyobject==L.pyobject, u.pyobject, bcs=bcs.pyobject)
export lvsolve

nlvsolve(F,u;bcs=nothing)=fenics.solve(F.pyobject==0,u.pyobject,bcs=bcs.pyobject)
export nlvsolve

anlvsolve(F,a,u,bcs,tol,M)=fenics.solve(F.pyobject==a.pyobject,u.pyobject,bcs=bcs.pyobject,tol=tol,M=M)
#this function hasnt been tested yet, so isnt exported
#TODO :test it

"""
errornorm is the function to calculate the error between our exact and calculated
solution. The norm kwarg defines the norm measure used (by default it is the L2 norm)
"""
errornorm(ans,sol;norm="L2") = fenics.errornorm(ans.pyobject,sol.pyobject,norm)
export errornorm
