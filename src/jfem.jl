
#These are the commands to define the Fem class, and assemble the Matrix in Julia
#Tests for these can be found in the TODO


using FEniCS


@fenicsclass FunctionSpace

FunctionSpace(mesh::Mesh, family::String, degree::Int) = FunctionSpace(fenics.FunctionSpace(mesh.pyobject, family, degree))
#add more functionspace functions

@fenicsclass Argument

TrialFunction(V::FunctionSpace) = Argument(fenics.TrialFunction(V.pyobject))
TestFunction(V::FunctionSpace) = Argument(fenics.TestFunction(V.pyobject))

@fenicsclass Constant

Constant(x::Real) = Constant(fenics.Constant(x, name="Constant($x)"))

@fenicsclass Expression
Expression(expr::String) = Expression(fenics.Expression(expr))
  
