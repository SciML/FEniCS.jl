
#These are the commands to define the Mesh class in Julia.
#Tests for these can be found in the test_create.jl and test_pycreate.jl

@fenicsclass Mesh  #https://fenicsproject.org/olddocs/dolfin/1.5.0/python/programmers-reference/cpp/mesh/Mesh.html
#possibly look at the type of object returned each time ( I believe numpy_arrays)
#are converted automatically by PyCall
cell_orientations(mesh::Mesh) = fenicspycall(mesh, :cell_orientations)
cells(mesh::Mesh) = fenicspycall(mesh, :cells)
hmin(mesh::Mesh) = fenicspycall(mesh, :hmin)
hmax(mesh::Mesh) = fenicspycall(mesh, :hmax)
init(mesh::Mesh) = fenicspycall(mesh, :init)
init(mesh::Mesh, dim::Int) = fenicspycall(mesh, :init, dim) # version with dims
init_global(mesh::Mesh) = fenicspycall(mesh,:init_global)
#mpi_comm
coordinates(mesh::Mesh) = fenicspycall(mesh,:coordinates)
#color
data(mesh::Mesh) = fenicspycall(mesh,:data)
domains(mesh::Mesh) = fenicspycall(mesh, :domains)
geometry(mesh::Mesh) = fenicspycall(mesh,:geometry)
num_cells(mesh::Mesh) = fenicspycall(mesh, :num_cells)
num_edges(mesh::Mesh) = fenicspycall(mesh, :num_edges)
num_entities(mesh::Mesh, dim::Int) = fenicspycall(mesh, :num_entities, dim)
num_faces(mesh::Mesh) = fenicspycall(mesh, :num_faces)
num_facets(mesh::Mesh) = fenicspycall(mesh, :num_facets)
num_vertices(mesh::Mesh) = fenicspycall(mesh, :num_vertices)
#hash(mesh::Mesh) = fenicspycall(Mesh, :hash)
bounding_box_tree(mesh::Mesh) = fenicspycall(mesh,:bounding_box_tree) #this object is a pyobject

export cell_orientations,cells,hmin , hmax, init, init_global, coordinates, data,
domains, geometry,num_cells,num_edges,num_entities,num_faces,num_facets,num_vertices, bounding_box_tree

UnitTriangleMesh() = Mesh(fenics.UnitTriangleMesh())

"""
  Mesh
Mesh is equivanlent to the Mesh function in fenics
"""

#name change
Mesh(path::String) = Mesh(fenics.Mesh(path))

UnitTetrahedronMesh() = Mesh(fenics.UnitTetrahedronMesh())

UnitSquareMesh(nx::Int, ny::Int, diagonal::String="right") = Mesh(fenics.UnitSquareMesh(nx, ny, diagonal))

UnitQuadMesh(nx::Int,ny::Int) = Mesh(fenics.UnitQuadMesh(nx,ny))

UnitIntervalMesh(nx::Int) = Mesh(fenics.UnitIntervalMesh(nx))

UnitCubeMesh(nx::Int, ny::Int, nz::Int) = Mesh(fenics.UnitCubeMesh(nx,ny,nz))

BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int)= Mesh(fenics.BoxMesh(p0,p1,nx,ny,nz)) #look at how to define fenics.point

RectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::String="right") = Mesh(fenics.RectangleMesh(p0,p1,nx,ny))

export UnitTriangleMesh, UnitTetrahedronMesh, UnitSquareMesh, UnitQuadMesh,
UnitIntervalMesh, UnitCubeMesh, BoxMesh, RectangleMesh, Mesh

function pyUnitTriangleMesh()
  pycall(fenics.UnitTriangleMesh::PyObject,PyObject::Type)
end

function pyUnitTetrahedronMesh()
  pycall(fenics.UnitTetrahedronMesh::PyObject,PyObject::Type)
end

function pyUnitCubeMesh(nx::Int, ny::Int, nz::Int)
  pycall(fenics.UnitCubeMesh::PyObject,PyObject::Type,nx,ny,nz)
end

function pyBoxMesh(p0, p1, nx::Int, ny::Int, nz::Int) # look at array types to declare p0,p1
  pycall(fenics.BoxMesh::PyObject,PyObject::Type,p0,p1,nx,ny,nz)
end

function pyRectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::String="right")
  pycall(fenics.RectangleMesh::PyObject,PyObject::Type,p0,p1,nx,ny,diagdir)
end

"""
For the diagdir, the possible options can be found below (these indicate the direction of the diagonals)
  (“left”, “right”, “right/left”, “left/right”, or “crossed”).
"""

function pyUnitSquareMesh(nx::Int,ny::Int,diagdir::String="right")
  pycall(fenics.UnitSquareMesh::PyObject,PyObject::Type,nx,ny,diagdir)
end

function pyUnitQuadMesh(nx::Int,ny::Int)
  pycall(fenics.UnitSquareMesh::PyObject,PyObject::Type,nx,ny)
end #https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/cpp/mesh/UnitQuadMesh.html
#states that the UnitQuadMesh code is experimental. Nevertheless I plan to add it , and maybe remove it at the final
#iteration

function pyUnitIntervalMesh(nx::Int)
  pycall(fenics.UnitIntervalMesh::PyObject,PyObject::Type,nx)
end

function pyMesh(path::String)
  pycall(fenics.Mesh::PyObject,PyObject::Type,path)
end

function Point(point::Vector) #a different data type was suggested. Will Investigate when I return to UK
  pycall(fenics.Point::PyObject,PyObject::Type,point)
end


export pyUnitTriangleMesh, pyUnitTetrahedronMesh, pyUnitSquareMesh, pyUnitQuadMesh,
pyUnitIntervalMesh, pyUnitCubeMesh, pyBoxMesh, pyRectangleMesh,pyMesh, Point
