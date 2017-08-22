# FEniCS

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/FEniCS.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/FEniCS.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/FEniCS.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiffEq/FEniCS.jl?branch=master)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/FEniCS.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/FEniCS.jl?branch=master)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations
of PDEs. This wrapper includes three parts:

1. Installation and direct access to FEniCS via a Conda installation. Alternatively one may use their current FEniCS installation https://fenicsproject.org/download/ .
2. A low-level development API and provides some functionality to make directly dealing with the library a little bit easier, but still requires knowledge of FEniCS itself. Interfaces have been provided for the main functions and their attributes, and instructions to add further ones have been provided here (https://gist.github.com/ysimillides/160bbf50a7e99d6656398aee48c88ef7)  .
3. A high-level API for usage with DifferentialEquations (https://github.com/JuliaDiffEq/DifferentialEquations.jl) which has not been implemented yet.

Various gists/jupyter notebooks show how to use the functionality(https://gist.github.com/ysimillides/d73c501f7f151c35d5002e8a11908a6e) .Any significant differences from FEniCS can be found here(https://gist.github.com/ysimillides/30fdbd71951c5e29e5f999635410b6de). It is used as an addon in the
DifferentialEquations.jl ecosystem. It is recommended that users use this
functionality through DifferentialEquations.jl. Paraview (https://www.paraview.org/) can also be used to visualize various results just like in FEniCS .

To get the wrapper on your system, follow the below steps:

1. Clone the package into Julia. This can be done via Pkg.clone("https://github.com/JuliaDiffEq/FEniCS.jl")
2. Proceed to build the package via Pkg.build("FEniCS")
3. Should be available to use. Due to the fact that FEniCS.jl uses the conda distribution to install and use FEniCS, this package will not currently work on Windows based systems due to the fact that a Windows Conda distribution is not currently supported by FEniCS ( https://fenicsproject.org/download/ ).

Note: Any suggestions/improvements/comments etc are always welcomed and can be made either on GitHub or via the gitter channel above.
This wrapper was originally started via the Google Summer of Code program (https://summerofcode.withgoogle.com/projects/#5988523772477440) along with the help of Chris Rackauckas and Bart Janssens.
