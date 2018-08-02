 #this file contains function related to plotting the mesh objects.
#we will use matplotlib.pyplot to achieve this
#
using FEniCS

import PyPlot: triplot, plot_trisurf, tripcolor, tricontourf

@pyimport matplotlib.tri as tri
@pyimport fenics

#creates a Triangulation for the Mesh type
function mesh2triangle(object::Mesh)
    xy = coordinates(object)
    return tri.Triangulation(xy[:,1],xy[:,2],cells(object))
end

#creates a Triangulation for the feMesh type
function mesh2triangle(space::feMesh)
  xy = space.nodes
  #needed to convert the element numbering back to a zero based index
  cells = space.elements.-1
  return tri.Triangulation(xy[:,1],xy[:,2],cells)
end



function plot(object::Mesh;kws...)
  geom = geometry(object)
  gdim = geom[:dim]()
  topol = topology(object)
  tdim = topol[:dim]()
  #plot 2D shapes
  if gdim == 2 && tdim == 2
    xy = coordinates(object)
    triangles = mesh2triangle(object)
    triplot(triangles; kws...)
  elseif gdim ==3 && tdim ==3
    bmesh = BoundaryMesh(object,"exterior",false)
    plot(bmesh;kws...)
  elseif gdim == 3 && tdim == 2
    xy = coordinates(object)
    plot_trisurf(xy[:,1],xy[:,2],xy[:,3],triangles = cells(object);kws...)
  end
end


#Only very basic function plotting currently exists, need to dispatch based on kwargs for full plotting.
function plot(object::Union{Expression,FeFunction};kws...)
  f_space = object.pyobject[:function_space]()
  mesh = f_space[:mesh]()
  geometry = mesh[:geometry]()
  topology = mesh[:topology]()
  gdim = geometry[:dim]()
  tdim = topology[:dim]()
  fvec = vector(object)
  if fvec[:size]() == mesh[:num_cells]()
    C = fvec[:array]()
    if gdim == 2 && tdim == 2
      tripcolor(mesh2triangle(mesh), C; kws...)
    end
  elseif gdim== 3 && tdim==2
    xy = coordinates(object)
    plot_trisurf(mesh2triangle(mesh),xy[:,2],C;kws...)
  elseif object.pyobject[:value_rank]()==0
    C = object.pyobject[:compute_vertex_values](mesh)
    tricontourf(mesh2triangle(Mesh(mesh)),C,40;kws...)
  end
end

function surf_plot(object::Union{Expression,FeFunction};kws...)
  f_space = object.pyobject[:function_space]()
  mesh = f_space[:mesh]()
  geometry = mesh[:geometry]()
  topology = mesh[:topology]()
  gdim = geometry[:dim]()
  tdim = topology[:dim]()
  fvec = vector(object)
  if object.pyobject[:value_rank]()==0
    C = object.pyobject[:compute_vertex_values](mesh)
    xy = coordinates(Mesh(mesh))
    PyPlot.surf(xy[:,1],xy[:,2],C)
  else
    error("surf_plot is not available for these dimensions")
  end
end


#plotting for solution computed via jinterface.jl, directly with the Mesh
function plot(mesh::Mesh, solution::AbstractArray, levels::Int=40;kws...)
  tricontourf(mesh2triangle(mesh),solution,levels;kws...)
end

#plotting for solution computed via jinterface.jl
function plot(space::feMesh,solution::AbstractArray, levels::Int=40;kws...)
  tricontourf(mesh2triangle(space),solution,levels;kws...)
end

export plot,surf_plot
