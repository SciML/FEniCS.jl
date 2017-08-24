using FEniCS
using Base.Test
using PyCall
#using PyPlot
#using Plots

@pyimport fenics

mesh = UnitSquareMesh(4,4)
V = FunctionSpace(mesh,"P",1)
u_D = Expression("1+x[0]*x[0]+2*x[1]*x[1]",degree=2 )
bc1 = DirichletBC(V,u_D, "on_boundary")
u=TrialFunction(V)
v=TestFunction(V)
m=Constant(-6.0)
a = inner(nabla_grad(u),nabla_grad(v))*dx
b = dot(grad(u),grad(v))*dx
A  = assemble(a)
B  = assemble(b)
L = m*v*dx
matrix1=get_array(A)
matrix2 =get_array(B)
U = FEniCS.Function(V)
lvsolve(a,L,U,bc1)
matrix3 = get_array(U)

#this testset checks that our wrapper to get the array returns an array of correct size.
@testset "get_array" begin
   @test matrix1 == matrix2
   @test size(matrix1) == (25,25)
   @test size(matrix3) == (25,)
 end


@test include("test_create.jl")
@test include("test_pycreate.jl")
@test include("test_jfem.jl")
#@test include("test_misc.jl")
