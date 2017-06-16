using FEniCS
using Base.Test

#ask about where the creation of be


test_cube = UnitCubeMesh(10,10,10)
x5 = hmin(test_cube)

# for roughly equals Julia uses \approx
#
# @testset "mesh" begin
#   @test ...
#   @test ...
#   @test...
# end
@test 1 != 2
@test 2 == 2
@test x5 ≈ 0.17320508075688767

@test include("test_create.jl")
@test include("test_pycreate.jl")
