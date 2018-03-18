
#this file includes the creation of the various julian meshes, returning true if all of
#them are creates succesfully with no errors

using FEniCS
using PyCall
@pyimport fenics


test_triangle = UnitTriangleMesh()
test_tetrahedron = UnitTetrahedronMesh()
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


try
  pyimport("mshr")
  circle1 = Circle(Point([0.0, 0.0]),3)
  circle2 = Circle(Point([0.0, 0.0]),4)
  rectangle = Rectangle(Point([0.0, 0.0]),Point([1.0, 1.0]))
  box = Box(Point([0.0, 0.0, 0.0]),Point([1.0, 1.0,1.0]))
  cone = Cone(Point([0.0,0.0,1.0]),Point([0.0, 0.0, 0.0]),3,50)
  sphere = Sphere(Point([0.0, 0.0, 0.0]),3)
  mesh1 = generate_mesh(circle1,64)
  mesh2= generate_mesh(circle2,64)
  circle3 = circle1+circle2
  circle4 = circle2-circle1
  mesh3 = generate_mesh(circle3,64)
  mesh4 = generate_mesh(circle4,64)
catch ee
 print("mshr has not been included, so it has not been tested")
end

true
