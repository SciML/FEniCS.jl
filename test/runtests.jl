using FEniCS
using Base.Test
using PyCall
using PyPlot
using Plots

@pyimport fenics


# for roughly equals Julia uses \approx
#
# @testset "mesh" begin
#   @test ...
#   @test ...
#   @test...
# end

@test include("test_create.jl")
@test include("test_pycreate.jl")
#@test include("test_fem.jl")
#@test include("test_misc.jl")
