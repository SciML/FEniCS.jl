using FEniCS
using Test

@testset "assign" begin
    mesh = UnitIntervalMesh(10)
    V = FunctionSpace(mesh, "P", 1)
    u = interpolate(Expression("x[0]", degree = 1), V)
    arr1 = get_array(u)
    arr2 = randn(size(arr1))
    assign(u, arr2)
    @test get_array(u) == arr2
    assign(u, arr1)
    @test get_array(u) == arr1
end
