
using FEniCS

mesh = UnitSquareMesh(8, 8)
element = FiniteElement("Lagrange", ufl_cell(mesh))

family_element = family(element)
cell_element = cell(element)
degree_element = degree(element)
variant_element = variant(element)
reconstructed_element = reconstruct(element)
sobolev_space_element = sobolev_space(element)
facet_element = FacetNormal(mesh)


V = FunctionSpace(mesh, "P", 1)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)

A_1 = inner(grad(u),grad(v))*dx
A_2 = dot(grad(u),grad(v))*dx
A_3 = cross(grad(u),grad(v))
A_4 = outer(grad(u),grad(v))
ng = nabla_grad(u)

B_1 = assemble(A_1)
B_2 = assemble(A_2)

measure_x = dx
measure_s = ds
measure_S = dS
measure_P = dP


true
