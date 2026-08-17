#type alias for string or symbol
StringOrSymbol = Union{String, Symbol}
"""
    Mesh

Abstract wrapper for a FEniCS/DOLFIN mesh.

Concrete values contain a `PyCall.PyObject` in the internal `pyobject` field.
Use the constructors below to load or generate meshes; extend mesh behavior on
`Mesh` rather than on its implementation type.
"""
abstract type Mesh <: fenicsobject end

struct MeshImpl <: Mesh
    pyobject::PyObject
end

Mesh(pyobject::PyObject) = MeshImpl(pyobject)
#are converted automatically by PyCall

"""
    cell_orientations(mesh::Mesh)

Return the orientation associated with each cell in `mesh`.

The result is supplied by the wrapped FEniCS mesh and is useful when a
finite-element computation needs the orientation data explicitly.

# Arguments

- `mesh`: Mesh whose cell orientations are requested.

# Returns

The FEniCS cell-orientation array.
"""
cell_orientations(mesh::Mesh) = fenicspycall(mesh, :cell_orientations)

"""
    cells(mesh::Mesh)

Return the cell-to-vertex connectivity of `mesh`.

# Arguments

- `mesh`: Mesh whose connectivity is requested.

# Returns

The FEniCS cell-connectivity array.
"""
cells(mesh::Mesh) = fenicspycall(mesh, :cells)

"""
    hmin(mesh::Mesh)

Return the minimum cell diameter in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.

# Returns

The minimum cell diameter as reported by FEniCS.
"""
hmin(mesh::Mesh) = fenicspycall(mesh, :hmin)

"""
    hmax(mesh::Mesh)

Return the maximum cell diameter in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.

# Returns

The maximum cell diameter as reported by FEniCS.
"""
hmax(mesh::Mesh) = fenicspycall(mesh, :hmax)

"""
    init(mesh::Mesh)
    init(mesh::Mesh, dim::Int)

Initialize mesh connectivity data.

The one-argument form initializes all connectivity data. The two-argument
form initializes connectivity involving topological dimension `dim`.

# Arguments

- `mesh`: Mesh whose connectivity should be initialized.
- `dim`: Optional topological dimension used to restrict initialization.
"""
init(mesh::Mesh) = fenicspycall(mesh, :init)
init(mesh::Mesh, dim::Int) = fenicspycall(mesh, :init, dim) # version with dims

"""
    init_global(mesh::Mesh)

Initialize global mesh connectivity data for `mesh`.

# Arguments

- `mesh`: Mesh whose global connectivity should be initialized.
"""
init_global(mesh::Mesh) = fenicspycall(mesh, :init_global)

"""
    coordinates(mesh::Mesh)

Return the coordinates of all vertices in `mesh`.

# Arguments

- `mesh`: Mesh whose vertex coordinates are requested.

# Returns

An array containing one coordinate vector per mesh vertex.
"""
coordinates(mesh::Mesh) = fenicspycall(mesh, :coordinates)

"""
    data(mesh::Mesh)

Return the auxiliary data object associated with `mesh`.

# Arguments

- `mesh`: Mesh whose auxiliary data is requested.
"""
data(mesh::Mesh) = fenicspycall(mesh, :data)

"""
    domains(mesh::Mesh)

Return the domain markers associated with `mesh`.

# Arguments

- `mesh`: Mesh whose domain markers are requested.
"""
domains(mesh::Mesh) = fenicspycall(mesh, :domains)

"""
    topology(mesh::Mesh)

Return the topological structure of `mesh`.

# Arguments

- `mesh`: Mesh whose topology object is requested.
"""
topology(mesh::Mesh) = fenicspycall(mesh, :topology)

"""
    geometry(mesh::Mesh)

Return the geometric structure of `mesh`.

# Arguments

- `mesh`: Mesh whose geometry object is requested.
"""
geometry(mesh::Mesh) = fenicspycall(mesh, :geometry)

"""
    num_cells(mesh::Mesh)

Return the number of cells in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
num_cells(mesh::Mesh) = fenicspycall(mesh, :num_cells)

"""
    num_edges(mesh::Mesh)

Return the number of edges in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
num_edges(mesh::Mesh) = fenicspycall(mesh, :num_edges)

"""
    num_entities(mesh::Mesh, dim::Int)

Return the number of mesh entities of topological dimension `dim`.

# Arguments

- `mesh`: Mesh to inspect.
- `dim`: Topological dimension of the entities to count.
"""
num_entities(mesh::Mesh, dim::Int) = fenicspycall(mesh, :num_entities, dim)

"""
    num_faces(mesh::Mesh)

Return the number of faces in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
num_faces(mesh::Mesh) = fenicspycall(mesh, :num_faces)

"""
    num_facets(mesh::Mesh)

Return the number of facets in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
num_facets(mesh::Mesh) = fenicspycall(mesh, :num_facets)

"""
    num_vertices(mesh::Mesh)

Return the number of vertices in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
num_vertices(mesh::Mesh) = fenicspycall(mesh, :num_vertices)

"""
    bounding_box_tree(mesh::Mesh)

Return the bounding-box tree associated with `mesh`.

# Arguments

- `mesh`: Mesh whose spatial index is requested.
"""
bounding_box_tree(mesh::Mesh) = fenicspycall(mesh, :bounding_box_tree) #this object is a pyobject

"""
    rmax(mesh::Mesh)

Return the maximum cell inradius in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
rmax(mesh::Mesh) = fenicspycall(mesh, :rmax)

"""
    rmin(mesh::Mesh)

Return the minimum cell inradius in `mesh`.

# Arguments

- `mesh`: Mesh to inspect.
"""
rmin(mesh::Mesh) = fenicspycall(mesh, :rmin)
"""
    size(mesh::Mesh, dim::Int)

Return the number of local mesh entities in topological dimension `dim`.

# Arguments
- `mesh`: FEniCS mesh to query.
- `dim`: Topological dimension of the requested entities.

# Examples
```julia
julia> size(mesh, 0) # vertices
4
```
"""
size(mesh::Mesh, dim::Int) = fenicspycall(mesh, :size, dim) # version with dims
"""
    ufl_cell(mesh::Mesh)

Return the UFL cell associated with `mesh`.
"""
ufl_cell(mesh::Mesh) = fenicspycall(mesh, :ufl_cell)

"""
    ufl_domain(mesh::Mesh)

Return the UFL domain associated with `mesh`.
"""
ufl_domain(mesh::Mesh) = fenicspycall(mesh, :ufl_domain)

"""
    ufl_id(mesh::Mesh)

Return the UFL identifier associated with `mesh`.
"""
ufl_id(mesh::Mesh) = fenicspycall(mesh, :ufl_id)

"""
    CellDiameter(mesh::Mesh)

Construct the symbolic cell-diameter expression for `mesh`.

# Arguments

- `mesh`: Mesh used to determine the cell diameter.

# Returns

A symbolic [`Expression`](@ref) suitable for use in a variational form.
"""
CellDiameter(mesh::Mesh) = Expression(fenics.CellDiameter(mesh.pyobject))

"""
    CellNormal(mesh::Mesh)

Construct the symbolic cell-normal expression for `mesh`.

# Arguments

- `mesh`: Mesh used to determine the cell normal.

# Returns

A symbolic [`Expression`](@ref) suitable for use in a variational form.
"""
CellNormal(mesh::Mesh) = Expression(fenics.CellNormal(mesh.pyobject))

"""
    CellVolume(mesh::Mesh)

Construct the symbolic cell-volume expression for `mesh`.

# Arguments

- `mesh`: Mesh used to determine the cell volume.

# Returns

A symbolic [`Expression`](@ref) suitable for use in a variational form.
"""
CellVolume(mesh::Mesh) = Expression(fenics.CellVolume(mesh.pyobject))

export cell_orientations, cells, hmin, hmax, init, init_global, coordinates, data,
    domains, geometry, topology, num_cells, num_edges, num_entities, num_faces,
    num_facets, num_vertices, bounding_box_tree,
    rmax, rmin, size, ufl_cell, ufl_domain, ufl_id, CellDiameter, CellNormal, CellVolume

# This constant is initialized in __init__
export CellType

"""
    Mesh(path::StringOrSymbol)

Load a FEniCS mesh from `path`.

# Arguments

- `path`: Filename or path understood by FEniCS.
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
UnitTriangleMesh() = Mesh(fenics.cpp.generation.UnitTriangleMesh.create())

"""
UnitTetrahedronMesh() \n
A mesh consisting of a single tetrahedron with vertices at \n
(0, 0, 0) (1, 0, 0) (0, 1, 0) (0, 0, 1)
"""
UnitTetrahedronMesh() = Mesh(fenics.cpp.generation.UnitTetrahedronMesh.create())

"""
UnitSquareMesh(nx::Int, ny::Int, diagonal::StringOrSymbol="right" ) \n

Triangular/quadrilateral mesh of the 2D unit square [0,1] x [0,1]. \n
Given the number of cells (nx, ny) in each direction, the total number of triangles \n
will be 2*nx*ny and the total number of vertices will be (nx + 1)*(ny + 1) \n
diagonal ("left", "right", "right//left", "left//right", or "crossed") indicates the direction of the diagonals.
"""
function UnitSquareMesh(nx::Int, ny::Int, diagonal::StringOrSymbol = "right")
    return Mesh(fenics.UnitSquareMesh(nx, ny, diagonal))
end
function UnitSquareMesh(nx::Int, ny::Int, cellType::PyObject)
    return Mesh(fenics.UnitSquareMesh.create(nx, ny, cellType))
end

"""
    UnitQuadMesh(nx::Int, ny::Int)

Deprecated compatibility helper for constructing a unit quadrilateral mesh.

Use a current FEniCS quadrilateral mesh constructor instead.
"""
function UnitQuadMesh(nx::Int, ny::Int)
    return println("Deprecated in FEniCS v.2018, remove in .7 Julia")
end

"""
    UnitIntervalMesh(nx::Int)

Construct a mesh of the unit interval `(0, 1)` with `nx` cells and `nx + 1`
vertices.

# Arguments

- `nx`: Number of cells in the interval.
"""
UnitIntervalMesh(nx::Int) = Mesh(fenics.UnitIntervalMesh(nx))

"""
    UnitCubeMesh(nx::Int, ny::Int, nz::Int)

Construct a tetrahedral mesh of the three-dimensional unit cube.

The mesh has `6 * nx * ny * nz` tetrahedra and
`(nx + 1) * (ny + 1) * (nz + 1)` vertices.

# Arguments

- `nx`, `ny`, `nz`: Number of cells in each coordinate direction.
"""
UnitCubeMesh(nx::Int, ny::Int, nz::Int) = Mesh(fenics.UnitCubeMesh(nx, ny, nz))
function UnitCubeMesh(nx::Int, ny::Int, nz::Int, cellType::PyObject)
    return Mesh(fenics.UnitCubeMesh.create(nx, ny, nz, cellType))
end
"""
    BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int)

Construct a tetrahedral mesh of the rectangular prism between `p0` and
`p1`.

# Arguments

- `p0`, `p1`: Opposite prism corners.
- `nx`, `ny`, `nz`: Number of cells in each coordinate direction.
"""
BoxMesh(p0, p1, nx::Int, ny::Int, nz::Int) = Mesh(fenics.BoxMesh(p0, p1, nx, ny, nz))
function BoxMesh(p::NTuple{2, PyObject}, n::NTuple{3, Int}, cellType::PyObject)
    return Mesh(fenics.BoxMesh.create(p, n, cellType))
end

"""
    RectangleMesh(p0, p1, nx::Int, ny::Int, diagdir::StringOrSymbol = "right")

Construct a triangular mesh of the rectangle between `p0` and `p1`.

# Arguments

- `p0`, `p1`: Opposite rectangle corners.
- `nx`, `ny`: Number of cells in each coordinate direction.
- `diagdir`: Diagonal orientation: `"left"`, `"right"`, `"right/left"`,
  `"left/right"`, or `"crossed"`.
"""
function RectangleMesh(p0, p1, nx::Int, ny::Int, diagdir::StringOrSymbol = "right")
    return Mesh(fenics.RectangleMesh(p0, p1, nx, ny, diagdir))
end
function RectangleMesh(p::NTuple{2, PyObject}, n::NTuple{2, Int}, cellType::PyObject)
    return Mesh(fenics.RectangleMesh.create(p, n, cellType))
end

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
function BoundaryMesh(mesh::Mesh, type_boundary::StringOrSymbol = "exterior", order = true)
    return Mesh(fenics.BoundaryMesh(mesh.pyobject, type_boundary, order))
end

export UnitTriangleMesh, UnitTetrahedronMesh, UnitSquareMesh, UnitQuadMesh,
    UnitIntervalMesh, UnitCubeMesh, BoxMesh, RectangleMesh, Mesh, BoundaryMesh

"""
    pyUnitTriangleMesh()

Return the underlying Python reference triangle mesh object.
"""
function pyUnitTriangleMesh()
    return fenics.cpp.generation.UnitTriangleMesh.create()
end

"""
    pyUnitTetrahedronMesh()

Return the underlying Python reference tetrahedron mesh object.
"""
function pyUnitTetrahedronMesh()
    return fenics.cpp.generation.UnitTetrahedronMesh.create()
end

"""
    pyUnitCubeMesh(nx::Int, ny::Int, nz::Int)

Construct the underlying Python unit-cube mesh object.

# Arguments

- `nx`, `ny`, `nz`: Number of cells in each coordinate direction.
"""
function pyUnitCubeMesh(nx::Int, ny::Int, nz::Int)
    return pycall(fenics.UnitCubeMesh::PyObject, PyObject::Type, nx, ny, nz)
end

"""
    pyBoxMesh(p0, p1, nx::Int, ny::Int, nz::Int)

Construct the underlying Python box mesh object.

# Arguments

- `p0`, `p1`: Opposite box corners.
- `nx`, `ny`, `nz`: Number of cells in each coordinate direction.
"""
function pyBoxMesh(p0, p1, nx::Int, ny::Int, nz::Int) # look at array types to declare p0,p1
    return pycall(fenics.BoxMesh::PyObject, PyObject::Type, p0, p1, nx, ny, nz)
end

"""
    pyRectangleMesh(p0, p1, nx::Int, ny::Int, diagdir::StringOrSymbol = "right")

Construct the underlying Python rectangle mesh object.

# Arguments

- `p0`, `p1`: Opposite rectangle corners.
- `nx`, `ny`: Number of cells in each coordinate direction.
- `diagdir`: Diagonal orientation passed to FEniCS.
"""
function pyRectangleMesh(p0, p1, nx::Int, ny::Int, diagdir::StringOrSymbol = "right")
    return pycall(fenics.RectangleMesh::PyObject, PyObject::Type, p0, p1, nx, ny, diagdir)
end

"""
    pyUnitSquareMesh(nx::Int, ny::Int, diagdir::StringOrSymbol = "right")

Construct the underlying Python unit-square mesh object. `diagdir` may be
`"left"`, `"right"`, `"right/left"`, `"left/right"`, or `"crossed"`.
"""
function pyUnitSquareMesh(nx::Int, ny::Int, diagdir::StringOrSymbol = "right")
    return pycall(fenics.UnitSquareMesh::PyObject, PyObject::Type, nx, ny, diagdir)
end

"""
    pyUnitQuadMesh(nx::Int, ny::Int)

Construct the underlying Python quadrilateral mesh object.
"""
function pyUnitQuadMesh(nx::Int, ny::Int)
    return pycall(fenics.UnitSquareMesh::PyObject, PyObject::Type, nx, ny)
end #https://fenicsproject.org/olddocs/dolfin/2016.2.0/python/programmers-reference/cpp/mesh/UnitQuadMesh.html
#states that the UnitQuadMesh code is experimental. Nevertheless I plan to add it , and maybe remove it at the final
#iteration

"""
    pyUnitIntervalMesh(nx::Int)

Construct the underlying Python unit-interval mesh object.

# Arguments

- `nx`: Number of cells in the interval.
"""
function pyUnitIntervalMesh(nx::Int)
    return pycall(fenics.UnitIntervalMesh::PyObject, PyObject::Type, nx)
end

"""
    pyMesh(path::StringOrSymbol)

Load the underlying Python FEniCS mesh object from `path`.
"""
function pyMesh(path::StringOrSymbol)
    return pycall(fenics.Mesh::PyObject, PyObject::Type, path)
end

"""
    Point(point::Union{Vector, Tuple})

Construct an underlying Python FEniCS point from a Julia vector or tuple.

# Arguments

- `point`: Coordinate vector or tuple.
"""
function Point(point::Union{Vector, Tuple})
    return pycall(fenics.Point::PyObject, PyObject::Type, point)
end

export pyUnitTriangleMesh, pyUnitTetrahedronMesh, pyUnitSquareMesh, pyUnitQuadMesh,
    pyUnitIntervalMesh, pyUnitCubeMesh, pyBoxMesh, pyRectangleMesh, pyMesh, Point

"""
    Cell

Abstract wrapper for a FEniCS mesh cell.
"""
abstract type Cell <: fenicsobject end

struct CellImpl <: Cell
    pyobject::PyObject
end

Cell(pyobject::PyObject) = CellImpl(pyobject)
"""
    Cell(mesh::MeshImpl, i::Int)

Wrap cell `i` from a FEniCS mesh.

# Arguments

- `mesh`: Underlying mesh implementation.
- `i`: Cell index understood by FEniCS.
"""
Cell(mesh::MeshImpl, i::Int) = Cell(fenics.Cell(mesh.pyobject, i))

"""
    get_vertex_coordinates(cell::Cell)

Return the coordinates of the vertices of `cell`.
"""
get_vertex_coordinates(cell::Cell) = fenicspycall(cell, :get_vertex_coordinates)

"""
    h(cell::Cell)

Return the greatest distance between two vertices of `cell`.
"""
h(cell::Cell) = fenicspycall(cell, :h)

"""
    midpoint(cell::Cell)

Return the midpoint of `cell`.
"""
midpoint(cell::Cell) = fenicspycall(cell, :midpoint)

"""
    volume(cell::Cell)

Return the volume of `cell`.
"""
volume(cell::Cell) = fenicspycall(cell, :volume)

export Cell, get_vertex_coordinates, h, midpoint, volume
