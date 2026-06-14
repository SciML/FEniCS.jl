using FEniCS
using Test

@test include("test_create.jl")
@test include("test_pycreate.jl")
@test include("test_jfem.jl")
