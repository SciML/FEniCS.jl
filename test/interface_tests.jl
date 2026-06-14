using FEniCS
using Test

#tests relating to interface.jl file
@testset "Interface" begin
    global mesh = UnitSquareMesh(2, 2)
    dbc(x, y) = 1
    global u = feMesh(mesh, dbc)
    @test (u.n_nodes == 9)
    @test (u.n_elements == 8)
    @test (u.n_internal_nodes == 1)
    @test (u.n_boundary_nodes == 8)
    @test (u.node_vals == [1, 1, 1, 1, 0, 1, 1, 1, 1])
end
