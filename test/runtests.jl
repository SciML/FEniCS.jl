using FEniCS
using Test

FEniCS.set_log_level(FEniCS.WARNING)

examples_dir = joinpath(@__DIR__, "..", "examples")
@assert ispath(examples_dir)
example_filenames = readdir(examples_dir)
@assert "tutorial1.jl" in example_filenames

@testset "Example $(filename)" for filename in example_filenames
    path = joinpath(examples_dir, filename)
    @test (include(path) ;true)
    # Test that PyCall is not used in the example.
    # Use of PyCall indicates, that we did not wrap enough
    # funcionality
    s = read(path, String)
    @test !occursin("PyCall", s)
end

@testset "Creation" begin
   @test include("test_create.jl")
   @test include("test_pycreate.jl")
   @test include("test_jfem.jl")
end;

@testset "assign" begin
    mesh = UnitIntervalMesh(10)
    V = FunctionSpace(mesh, "P", 1)
    u = interpolate(Expression("x[0]", degree=1), V)
    arr1 = get_array(u)
    arr2 = randn(size(arr1))
    assign(u, arr2)
    @test get_array(u) == arr2
    assign(u, arr1)
    @test get_array(u) == arr1
end

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
