# FEniCS.jl: Finite Element PDE Solving in Julia

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/FEniCS/stable/)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations
of PDEs. This wrapper includes three parts:

1. Installation and direct access to FEniCS via a Conda installation. Alternatively, one may use their current [FEniCS installation](https://fenicsproject.org/download/).
2. A low-level development API and provides some functionality to make directly dealing with the library a little bit easier, but still requires knowledge of FEniCS itself. Interfaces have been provided for the main functions and their attributes, and instructions to add further ones can be [found here](https://gist.github.com/ysimillides/160bbf50a7e99d6656398aee48c88ef7).
3. A high-level API for usage with [DifferentialEquations](https://github.com/JuliaDiffEq/DifferentialEquations.jl). An example can be seen [solving the heat equation with high order adaptive timestepping](https://github.com/JuliaDiffEq/FEniCS.jl/blob/master/examples/heat_equation.jl).

Various [gists/jupyter](https://gist.github.com/ysimillides/20761c511a8807ae11c2b9e70606985e) notebooks have been created to provide a brief overview of the overall functionality, and of any differences between the Pythonic FEniCS and the Julian wrapper.
DifferentialEquations.jl ecosystem. [Paraview](https://www.paraview.org/) can also be used to visualize various results, just like in FEniCS (see below).

## Installation Instructions

To get the wrapper on your system, providing a FEniCS installation exists, follow the below steps:

1. Add PyCall with the correct Python environment corresponding to FEniCS. Then simply add FEniCS.jl using Pkg.add("FEniCS")

2. Alternatively, one can install [Docker](https://www.docker.com/) and then run the following command

```sh
docker run -ti cmhyett/julia-fenics
```
and once inside, Julia can be accessed by calling
```sh
julia
```
Once inside the Julia environment, simply add FEniCS with Pkg.add("FEniCS"). All other dependencies are handled by the docker image.

This wrapper was originally started via the [Google Summer of Code program](https://summerofcode.withgoogle.com/projects/#5988523772477440) along with the help of Chris Rackauckas and Bart Janssens. This was continued via [GSoC '18](https://summerofcode.withgoogle.com/projects/#6466456292687872) along with the help of Chris Rackauckas and Timo Betcke.

## Tutorial

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
U = FeFunction(V)
lvsolve(a,L,U,bc1) #linear variational solver
errornorm(u_D, U, norm="L2")
get_array(L) #this returns an array for the stiffness matrix
get_array(U) #this returns an array for the solution values
vtkfile = File("poisson/solution.pvd")
vtkfile << U.pyobject #exports the solution to a vtkfile
```

We can also plot the solution (this relies on FEniCS backend for plotting) or import it from our file into Paraview:

```julia
import PyPlot # plotting won't work if PyPlot is not imported
FEniCS.Plot(U)
FEniCS.Plot(mesh)

```

![alt text](https://user-images.githubusercontent.com/16087601/34915339-b77e8694-f91c-11e7-9ae1-db1e114a177a.png "Solution")

![alt text](https://user-images.githubusercontent.com/16087601/34915337-b2c0aede-f91c-11e7-986a-5658d23c262e.png "Square Mesh")

See the examples directory for more examples.

