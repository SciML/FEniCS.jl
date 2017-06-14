@fenicsclass Mesh  #https://fenicsproject.org/olddocs/dolfin/1.5.0/python/programmers-reference/cpp/mesh/Mesh.html
hmin(mesh::Mesh) = fenicspycall(mesh, :hmin)
hmax(mesh::Mesh) = fenicspycall(mesh, :hmax)
coordinates(mesh::Mesh) = fenicspycall(mesh,:coordinates)
bounding_box_tree(mesh::Mesh) = fenicspycall(mesh,:bounding_box_tree) #this object is a pyobject
export hmin,hmax,coordinates,bounding_box_tree

UnitTriangleMesh() = Mesh(fenics.UnitTriangleMesh())

UnitTetrahedronMesh() = Mesh(fenics.UnitTetrahedronMesh())

UnitSquareMesh(nx::Int, ny::Int, diagonal::String="right") = Mesh(fenics.UnitSquareMesh(nx, ny, diagonal))

UnitQuadMesh(nx::Int,ny::Int) = Mesh(fenics.UnitSquareMesh(nx,ny))

UnitIntervalMesh(nx::Int) = Mesh(fenics.UnitIntervalMesh(nx))

UnitCubeMesh(nx::Int, ny::Int, nz::Int) = Mesh(fenics.UnitCubeMesh(nx,ny,nz))

BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int)= Mesh(fenics.BoxMesh(p0,p1,nx,ny,nz)) #look at how to define fenics.point

RectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::String="right") = Mesh(fenics.RectangleMesh(p0,p1,nx,ny))


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
