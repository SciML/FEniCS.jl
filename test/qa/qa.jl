using FEniCS, Aqua, JET, Test

@testset "Aqua" begin
    Aqua.test_all(FEniCS)
end

@testset "JET" begin
    JET.test_package(FEniCS; target_defined_modules = true)
end
