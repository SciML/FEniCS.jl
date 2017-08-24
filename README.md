# FEniCS

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/FEniCS.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/FEniCS.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/FEniCS.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiffEq/FEniCS.jl?branch=master)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/FEniCS.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/FEniCS.jl?branch=master)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations
of PDEs. This wrapper includes three parts:

1. Installation and direct access to FEniCS via a Conda installation. Alternatively one may use their current [FEniCS installation](https://fenicsproject.org/download/) .
2. A low-level development API and provides some functionality to make directly dealing with the library a little bit easier, but still requires knowledge of FEniCS itself. Interfaces have been provided for the main functions and their attributes, and instructions to add further ones can be [found here](https://gist.github.com/ysimillides/160bbf50a7e99d6656398aee48c88ef7).
3. A high-level API for usage with [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl) which has not been implemented yet.

Various [gists/jupyter](https://gist.github.com/ysimillides/20761c511a8807ae11c2b9e70606985e) notebooks have been created to provide a brief overview of the overall functionality, and of any differences between the pythonic FEniCS and the julian wrapper.
DifferentialEquations.jl ecosystem. It is recommended that users use this
functionality through DifferentialEquations.jl.
[Paraview](https://www.paraview.org/) can also be used to visualize various results just like in FEniCS (see below) .


To get the wrapper on your system, follow the below steps:

1. Clone the package into Julia. This can be done via Pkg.clone("https://github.com/JuliaDiffEq/FEniCS.jl") and when registered via Pkg.add("FEniCS")
2. Proceed to build the package via Pkg.build("FEniCS")
3. Should be available to use. Due to the fact that FEniCS.jl uses the conda distribution to install and use FEniCS, the installation provided will not currently work on Windows based systems due to the fact that a Windows Conda distribution is not currently supported by [FEniCS](https://fenicsproject.org/download/).

Note: Any suggestions/improvements/comments etc are always welcomed and can be made either on GitHub or via the gitter channel above.
This wrapper was originally started via the [Google Summer of Code program](https://summerofcode.withgoogle.com/projects/#5988523772477440) along with the help of Chris Rackauckas and Bart Janssens.


Below is a small demonstration of how a user would use our code to solve the Poisson equation with Dirichlet conditions. This directly mirrors one of the **[tutorials](https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft01_poisson.py)** FEniCS provides 
```julia
using FEniCS
mesh = UnitSquareMesh(8,8) 
V = FunctionSpace(mesh,"P",1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]", degree=2)
u = TrialFunction(V)
bc1 = DirichletBC(V,u_D, "on_boundary")
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u),grad(v))*dx
L = f*v*dx
U = FEniCS.Function(V)
lvsolve(a,L,U,bc1) #linear variational solver
errornorm(u_D, U, norm="L2")
get_array(L) #this returns an array for the stiffness matrix
get_array(U) #this returns an array for the solution values
vtkfile = File('poisson/solution.pvd')
vtkfile << U.pyobject #exports the solution to a vtkfile
```

We can also plot the solution (this relies on FEniCS backend for plotting) or import it from our file into Paraview:

```julia
FEniCS.Plot(U)
FEniCS.Plot(mesh)

```

![alt text](https://github.com/ysimillides/FEniCS.jl/blob/master/examples/result.png "Solution")
 
![alt text](https://github.com/ysimillides/FEniCS.jl/blob/master/examples/mesh.png "Square Mesh")
