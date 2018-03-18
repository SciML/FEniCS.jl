using FEniCS
using Base.Test


@testset "Tutorials" begin
   @test include("tutorial1.jl")
   @test include("tutorial2.jl")
   @test include("tutorial3.jl")
   @test include("tutorial6.jl")
end;

@testset "Creation" begin
   @test include("test_create.jl")
   @test include("test_pycreate.jl")
   @test include("test_jfem.jl")
end;

print("tests completed")
