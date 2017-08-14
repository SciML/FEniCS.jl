
#this file includes the creation of the various julian meshes, returning true if all of
#them are creates succesfully with no errors

using FEniCS
using PyCall
@pyimport fenics


test_triangle = UnitTriangleMesh()
test_tetrahedron = UnitTetrahedronMesh()
#TODO test_quadmesh =UnitQuadMesh(10,10) # is one of the experimental meshes, and doesnt
#current work. I believe it needs an MPI_Comm?
test_interval = UnitIntervalMesh(10)
test_square = UnitSquareMesh(10,10,"crossed")
test_cube = UnitCubeMesh(10,10,10)
test_box = BoxMesh(Point([0.0,0.0,0.0]),Point([1.0,1.0,1.0]),10,10,10)
test_rectangle = RectangleMesh(Point([0.0, 0.0]), Point([10.0, 4.0]), 10, 10)

#the below functions simply check the creation of the objects,
#without (currently) verifying values
x1 = cell_orientations(test_triangle)
x2 = cells(test_triangle)
x3 = hmin(test_triangle)
x4 = hmax(test_triangle)
x5 = coordinates(test_triangle)
x6 = data(test_triangle)
x7 = domains(test_triangle)
x8 = geometry(test_triangle)
x9 = num_cells(test_triangle)
x10 = num_edges(test_triangle)
x11 = num_entities(test_triangle,1)
x12 = num_faces(test_triangle)
x13 = num_facets(test_triangle)
x14 = num_vertices(test_triangle)
x15 = bounding_box_tree(test_triangle)
x16 = init(test_triangle)
x17 = init(test_triangle,1)
#TODO : create a test file that checks that the various functions provide the correct results
true
#coordinates.jl (for the point)
#rotations.jl
