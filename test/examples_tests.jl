using FEniCS
using Test

FEniCS.set_log_level(FEniCS.WARNING)

examples_dir = joinpath(@__DIR__, "..", "examples")
@assert ispath(examples_dir)
example_filenames = readdir(examples_dir)
@assert "tutorial1.jl" in example_filenames

@testset "Example $(filename)" for filename in example_filenames
    path = joinpath(examples_dir, filename)
    @test (include(path); true)
    # Test that PyCall is not used in the example.
    # Use of PyCall indicates, that we did not wrap enough
    # functionality
    s = read(path, String)
    @test !occursin("PyCall", s)
end
