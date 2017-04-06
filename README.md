# FEniCS

[![Build Status](https://travis-ci.org/ChrisRackauckas/FEniCS.jl.svg?branch=master)](https://travis-ci.org/ChrisRackauckas/FEniCS.jl)

[![Coverage Status](https://coveralls.io/repos/ChrisRackauckas/FEniCS.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ChrisRackauckas/FEniCS.jl?branch=master)

[![codecov.io](http://codecov.io/github/ChrisRackauckas/FEniCS.jl/coverage.svg?branch=master)](http://codecov.io/github/ChrisRackauckas/FEniCS.jl?branch=master)

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
3. Should be available to use. Due to errors with PyCall, it does not currently work on a Windows based system.

Note: This is a work in progress!
