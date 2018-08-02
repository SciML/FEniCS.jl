#type alias for string or symbol
StringOrSymbol = Union{String,Symbol}
@fenicsclass Mesh  #https://fenicsproject.org/olddocs/dolfin/1.5.0/python/programmers-reference/cpp/mesh/Mesh.html
#are converted automatically by PyCall

# Get the cell orientations set.
cell_orientations(mesh::Mesh) = fenicspycall(mesh, :cell_orientations)
#returns cell connectivity
cells(mesh::Mesh) = fenicspycall(mesh, :cells)
#Compute minimum cell diameter.
hmin(mesh::Mesh) = fenicspycall(mesh, :hmin)
#Compute maximum cell diameter.
hmax(mesh::Mesh) = fenicspycall(mesh, :hmax)
init(mesh::Mesh) = fenicspycall(mesh, :init)
init(mesh::Mesh, dim::Int) = fenicspycall(mesh, :init, dim) # version with dims
init_global(mesh::Mesh) = fenicspycall(mesh,:init_global)
#returns coordinates of all vertices
coordinates(mesh::Mesh) = fenicspycall(mesh,:coordinates)

data(mesh::Mesh) = fenicspycall(mesh,:data)
#Get  mesh (sub)domains
domains(mesh::Mesh) = fenicspycall(mesh, :domains)
topology(mesh::Mesh) = fenicspycall(mesh,:topology)

#get mesh geometry
geometry(mesh::Mesh) = fenicspycall(mesh,:geometry)
#returns number of cells
num_cells(mesh::Mesh) = fenicspycall(mesh, :num_cells)
#returns number of edges
num_edges(mesh::Mesh) = fenicspycall(mesh, :num_edges)
#Get number of entities of given topological dimension.
num_entities(mesh::Mesh, dim::Int) = fenicspycall(mesh, :num_entities, dim)
#Get number of faces in mesh.
num_faces(mesh::Mesh) = fenicspycall(mesh, :num_faces)
#Get number of facets in mesh.
num_facets(mesh::Mesh) = fenicspycall(mesh, :num_facets)
#Get number of vertices in mesh.
num_vertices(mesh::Mesh) = fenicspycall(mesh, :num_vertices)
#hash(mesh::Mesh) = fenicspycall(Mesh, :hash)
bounding_box_tree(mesh::Mesh) = fenicspycall(mesh,:bounding_box_tree) #this object is a pyobject
#Compute maximum cell inradius.
rmax(mesh::Mesh)=fenicspycall(mesh, :rmax)
#Compute minimum cell inradius.
rmin(mesh::Mesh)=fenicspycall(mesh, :rmin)
#Get number of local entities of given topological dimension.
size(mesh::Mesh, dim::Int) = fenicspycall(mesh, :size, dim) # version with dims
#Returns the ufl cell of the mesh.
ufl_cell(mesh::Mesh)=fenicspycall(mesh, :ufl_cell)
#Returns the ufl Domain corresponding to the mesh.
ufl_domain(mesh::Mesh)=fenicspycall(mesh, :ufl_domain)
#Returns an id that UFL can use to decide if two objects are the same.
ufl_id(mesh::Mesh)=fenicspycall(mesh, :ufl_id)

export cell_orientations,cells,hmin , hmax, init, init_global, coordinates, data,
domains, geometry,topology, num_cells,num_edges,num_entities,num_faces,num_facets,num_vertices, bounding_box_tree,
rmax, rmin, size, ufl_cell , ufl_domain, ufl_id

"""
Mesh(path::StringOrSymbol) \n
Creates a Mesh based on a specified filename(path)
"""
Mesh(path::StringOrSymbol) = Mesh(fenics.Mesh(path))

"""
Mesh(object::Mesh) \n
Creates a copy of a mesh
"""
Mesh(object::Mesh) = Mesh(fenics.Mesh(object.pyobject))

"""
UnitTriangleMesh() \n
A mesh consisting of a single triangle with vertices at \n
(0, 0) (1, 0) (0, 1)
"""
UnitTriangleMesh() = Mesh(fenics.UnitTriangleMesh())

"""
UnitTetrahedronMesh() \n
A mesh consisting of a single tetrahedron with vertices at \n
(0, 0, 0) (1, 0, 0) (0, 1, 0) (0, 0, 1)
"""
UnitTetrahedronMesh() = Mesh(fenics.UnitTetrahedronMesh())


"""
UnitSquareMesh(nx::Int, ny::Int, diagonal::StringOrSymbol="right" ) \n

Triangular/quadrilateral mesh of the 2D unit square [0,1] x [0,1]. \n
Given the number of cells (nx, ny) in each direction, the total number of triangles \n
will be 2*nx*ny and the total number of vertices will be (nx + 1)*(ny + 1) \n
diagonal ("left", "right", "right//left", "left//right", or "crossed") indicates the direction of the diagonals.
"""
UnitSquareMesh(nx::Int, ny::Int, diagonal::StringOrSymbol="right") = Mesh(fenics.UnitSquareMesh(nx, ny, diagonal))


function UnitQuadMesh(nx::Int,ny::Int)
    println("Deprecated in FEniCS v.2018, remove in .7 Julia")
end

"""
UnitIntervalMesh(nx::Int) \n

A mesh of the unit interval (0, 1) with a given number of cells (nx) in the axial direction. \n
The total number of intervals will be nx and the total number of vertices will be (nx + 1).
"""
UnitIntervalMesh(nx::Int) = Mesh(fenics.UnitIntervalMesh(nx))


"""
UnitCubeMesh(nx::Int, ny::Int, nz::Int) \n

Tetrahedral/hexahedral mesh of the 3D unit cube [0,1] x [0,1] x [0,1]. \n
Given the number of cells (nx, ny, nz) in each direction, the total number of \n
tetrahedra will be 6*nx*ny*nz and the total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).

"""
UnitCubeMesh(nx::Int, ny::Int, nz::Int) = Mesh(fenics.UnitCubeMesh(nx,ny,nz))
"""
BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int) \n

Tetrahedral mesh of the 3D rectangular prism spanned by two points p0 and p1. \n
Given the number of cells (nx, ny, nz) in each direction, the total number of \n
tetrahedra will be 6*nx*ny*nz and the total number of vertices will be (nx + 1)*(ny + 1)*(nz + 1).
"""
BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int)= Mesh(fenics.BoxMesh(p0,p1,nx,ny,nz))

"""
RectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::StringOrSymbol="right") \n
Triangular mesh of the 2D rectangle spanned by two points p0 and p1. \n
Given the number of cells (nx, ny) in each direction, the total number \n
of triangles will be 2*nx*ny and the total number of vertices will be (nx + 1)*(ny + 1) \n
diagdir ("left", "right", "right/left", "left/right", or "crossed") indicates the direction of the diagonals.
"""
RectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::StringOrSymbol="right") = Mesh(fenics.RectangleMesh(p0,p1,nx,ny,diagdir))

"""
BoundaryMesh(mesh::Mesh,type_boundary::StringOrSymbol="exterior",order=true) \n

A BoundaryMesh  is a mesh over the boundary of some given mesh. \n
The cells of the boundary mesh (facets of the original mesh) are oriented to \n
produce outward pointing normals relative to the original mesh. \n
The type_boundary can be "exterior", "interior" or "local". "exterior" is the globally \n
external boundary, "interior" is the inter-process mesh and "local" is the boundary \n
of the local (this process) mesh. \n
order:(bool) Optional argument which can be used to control whether or not the \n
boundary mesh should be ordered according to the UFC ordering convention. \n
If set to false, the boundary mesh will be ordered with right-oriented facets \n
(outward-pointing unit normals). The default value is true.
"""
BoundaryMesh(mesh::Mesh,type_boundary::StringOrSymbol="exterior",order=true) = Mesh(fenics.BoundaryMesh(mesh.pyobject, type_boundary, order))

export UnitTriangleMesh, UnitTetrahedronMesh, UnitSquareMesh, UnitQuadMesh,
UnitIntervalMesh, UnitCubeMesh, BoxMesh, RectangleMesh, Mesh, BoundaryMesh

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

function pyRectangleMesh(p0,p1,nx::Int,ny::Int,diagdir::StringOrSymbol="right")
  pycall(fenics.RectangleMesh::PyObject,PyObject::Type,p0,p1,nx,ny,diagdir)
end

"""
For the diagdir, the possible options can be found below (these indicate the direction of the diagonals)
  (“left”, “right”, “right/left”, “left/right”, or “crossed”).
"""
function pyUnitSquareMesh(nx::Int,ny::Int,diagdir::StringOrSymbol="right")
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

function pyMesh(path::StringOrSymbol)
  pycall(fenics.Mesh::PyObject,PyObject::Type,path)
end

function Point(point::Union{Vector,Tuple})
  pycall(fenics.Point::PyObject,PyObject::Type,point)
end


export pyUnitTriangleMesh, pyUnitTetrahedronMesh, pyUnitSquareMesh, pyUnitQuadMesh,
pyUnitIntervalMesh, pyUnitCubeMesh, pyBoxMesh, pyRectangleMesh,pyMesh, Point
