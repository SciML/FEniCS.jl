using SafeTestsets

@safetestset "Explicit Imports" begin
    include("explicit_imports.jl")
end

@safetestset "Examples" begin
    include("examples_tests.jl")
end

@safetestset "Creation" begin
    include("creation_tests.jl")
end

@safetestset "assign" begin
    include("assign_tests.jl")
end

@safetestset "Interface" begin
    include("interface_tests.jl")
end
