# FEniCS

[![Join the chat at https://gitter.im/JuliaDiffEq/Lobby](https://badges.gitter.im/JuliaDiffEq/Lobby.svg)](https://gitter.im/JuliaDiffEq/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://travis-ci.org/JuliaDiffEq/FEniCS.jl.svg?branch=master)](https://travis-ci.org/JuliaDiffEq/FEniCS.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaDiffEq/FEniCS.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaDiffEq/FEniCS.jl?branch=master)
[![codecov.io](http://codecov.io/github/ChrisRackauckas/FEniCS.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/FEniCS.jl?branch=master)

#### Uncomment these if/when this is registered
#[![FEniCS](http://pkg.julialang.org/badges/FEniCS_0.5.svg)](http://pkg.julialang.org/?pkg=FEniCS)
#[![FEniCS](http://pkg.julialang.org/badges/FEniCS_0.6.svg)](http://pkg.julialang.org/?pkg=FEniCS)

FEniCS.jl is a wrapper for the FEniCS library for finite element discretizations
of PDEs. This wrapper includes three parts:

1. Installation and direct access to FEniCS
2. A low-level development API and provides some functionality to make directly dealing with the library a little bit easier, but still requires knowledge of FEniCS itself.
3. A high-level API for usage with DifferentialEquations

The tests show how to use the functionality. It is used as an addon in the
DifferentialEquations.jl ecosystem. It is recommended that users use this
functionality through DifferentialEquations.jl.

To get the wrapper on your system, follow the below steps:

1. Clone the package into Julia. This can be done via Pkg.clone("https://github.com/JuliaDiffEq/FEniCS.jl")
2. Proceed to build the package via Pkg.build("FEniCS")
3. Should be available to use. Due to the fact that FEniCS.jl uses the conda distribution to install and use FEniCS, this package will not currently work on Windows based systems due to the fact that a Windows Conda distribution is not currently supported by FEniCS ( https://fenicsproject.org/download/ )

Note: This is a work in progress, and the structure may change throughout the summer. Nevertheless any suggestions/improvements/comments etc are always welcomed.
