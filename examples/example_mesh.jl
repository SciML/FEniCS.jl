using FEniCS

"""
to create our mesh in both cases, we need to call the Julia function, and if necessary
(as in the case of the UnitCubeMesh provide it with some args. This is the same way we
would create it directly in FEniCS. The mesh returned is a Julian Object.
Should we wish to have it returned as a PyObject, functionality for that also exists ,
with the function name prefixed with py, eg pyUnitTriangleMesh.)
"""
test_triangle = UnitTriangleMesh()
test_cube = UnitCubeMesh(10,10,10)

test_pytriangle = pyUnitTriangleMesh()
test_pycube = pyUnitCubeMesh(10,10,10)

"""
Some of our functions also have \*optional\* args, such as the UnitSquareMesh.
By default, the diagdir (diagonal direction) is set as right, but may be adjusted
to other directions aswell (“left”, “right”, “right/left”, “left/right”, or “crossed”).
"""

test_square = UnitSquareMesh(10,10) # this will have a default direction of right
test_square2 = UnitSquareMesh(10,10,"left") # this has a direction of left.
"""
likewise, the PyObject meshes are called in the same way
"""
test_pysquare = pyUnitSquareMesh(10,10) # this will have a default direction of right
test_pysquare2 = pyUnitSquareMesh(10,10,"left") # this has a direction of left.

"""
To call various attributes (is this the correct word?) of the mesh class
depending on the type of the object (julian / PyObject) we call them in
slightly different ways. This is due to the nature of PyCall. For the julian object
we call it as foo(mesh), whereas for the pyobject we call it as mesh[:foo]().
Examples can be seen below
"""

x_1 = hmin(test_square)
x_2 = test_pysquare[:hmin]()
print(x_2)
#as a sanity check, we can always check that x_1 == x_2
print(x_1 == x_2)
