using FEniCS
using Test


@testset "Tutorials" begin
   @test include("tutorial1.jl")
   @test include("tutorial2.jl")
   @test include("tutorial3.jl")
   @test include("tutorial4.jl")
   @test include("tutorial5.jl")
   @test include("tutorial6.jl")
   @test include("tutorial7.jl")
   @test include("tutorial8.jl")
   @test include("tutorial9.jl")
   @test include("acoustic.jl")
   @test include("acoustic_new_assemble.jl")
end;

@testset "Creation" begin
   @test include("test_create.jl")
   @test include("test_pycreate.jl")
   @test include("test_jfem.jl")
end;

#tests relating to interface.jl file
@testset "Interface" begin
   global mesh = UnitSquareMesh(2,2)
   dbc(x,y) = 1
   global u = feMesh(mesh,dbc)
   @test (u.n_nodes==9)
   @test (u.n_elements==8)
   @test (u.n_internal_nodes==1)
   @test (u.n_boundary_nodes==8)
   @test (u.node_vals == [1,1,1,1,0,1,1,1,1])
end;


print("tests completed")
