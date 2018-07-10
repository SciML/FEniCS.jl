#This .jl files contains test relating to the FiniteElement class found in jfem.jl
#it currently only naively tests the creation of these.
#TODO write unit test for this
x=FiniteElement("Lagrange",fenics.triangle,1)
y=sobolev_space(x)
z=variant(x)
print(y)
print(z)
yy = family(x)
print("cell")
print(yy)
