using ExplicitImports
using FEniCS
using Test

@testset "ExplicitImports" begin
    @test check_no_implicit_imports(FEniCS) === nothing
    @test check_no_stale_explicit_imports(FEniCS) === nothing
end
