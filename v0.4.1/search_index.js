var documenterSearchIndex = {"docs":
[{"location":"#FEniCS.jl:-Finite-Element-PDE-Solving-in-Julia","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"","category":"section"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"(Image: Join the chat at https://julialang.zulipchat.com #sciml-bridged) (Image: Global Docs)","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"(Image: ColPrac: Contributor's Guide on Collaborative Practices for Community Packages) (Image: SciML Code Style)","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations of PDEs. This wrapper includes three parts:","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Installation and direct access to FEniCS via a Conda installation. Alternatively one may use their current FEniCS installation.\nA low-level development API and provides some functionality to make directly dealing with the library a little bit easier, but still requires knowledge of FEniCS itself. Interfaces have been provided for the main functions and their attributes, and instructions to add further ones can be found here.\nA high-level API for usage with DifferentialEquations. An example can be seen solving the heat equation with high order adaptive timestepping.","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Various gists/jupyter notebooks have been created to provide a brief overview of the overall functionality, and of any differences between the pythonic FEniCS and the julian wrapper. DifferentialEquations.jl ecosystem. Paraview can also be used to visualize various results just like in FEniCS (see below).","category":"page"},{"location":"#Installation-Instructions","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"Installation Instructions","text":"","category":"section"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"To get the wrapper on your system,providing a FEniCS installation exists, follow the below steps:","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Add PyCall with the correct python environment corresponding to FEniCS. Then simply add FEniCS.jl using Pkg.add(\"FEniCS\")\nAlternatively, one can install Docker and then run the following command","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"docker run -ti ysimillides/fenics-julia-docker","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"and once inside, 'julia' can be accessed by calling","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"julia","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Once inside the julia environment, simply add FEniCS with Pkg.add(\"FEniCS\"). All other dependencies are handled by the docker image.","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Note: Any suggestions/improvements/comments etc are always welcomed and can be made either on GitHub or via the gitter channel above. This wrapper was originally started via the Google Summer of Code program along with the help of Chris Rackauckas and Bart Janssens. This was continued via GSoC '18 along with the help of Chris Rackauckas and Timo Betcke.","category":"page"},{"location":"#Tutorial","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"Tutorial","text":"","category":"section"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"Below is a small demonstration of how a user would use our code to solve the Poisson equation with Dirichlet conditions. This directly mirrors one of the tutorials FEniCS provides","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"using FEniCS\nmesh = UnitSquareMesh(8,8)\nV = FunctionSpace(mesh,\"P\",1)\nu_D = Expression(\"1+x[0]*x[0]+2*x[1]*x[1]\", degree=2)\nu = TrialFunction(V)\nbc1 = DirichletBC(V,u_D, \"on_boundary\")\nv = TestFunction(V)\nf = Constant(-6.0)\na = dot(grad(u),grad(v))*dx\nL = f*v*dx\nU = FeFunction(V)\nlvsolve(a,L,U,bc1) #linear variational solver\nerrornorm(u_D, U, norm=\"L2\")\nget_array(L) #this returns an array for the stiffness matrix\nget_array(U) #this returns an array for the solution values\nvtkfile = File(\"poisson/solution.pvd\")\nvtkfile << U.pyobject #exports the solution to a vtkfile","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"We can also plot the solution (this relies on FEniCS backend for plotting) or import it from our file into Paraview:","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"import PyPlot # plotting won't work if PyPlot is not imported\nFEniCS.Plot(U)\nFEniCS.Plot(mesh)\n","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"(Image: alt text)","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"(Image: alt text)","category":"page"},{"location":"","page":"FEniCS.jl: Finite Element PDE Solving in Julia","title":"FEniCS.jl: Finite Element PDE Solving in Julia","text":"See the examples directory for more examples.","category":"page"}]
}