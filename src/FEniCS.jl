__precompile__(false)
module FEniCS
using PyCall
using Requires

fenics = pyimport_conda("fenics", "fenics=2019.1.0", "conda-forge")
ufl = pyimport_conda("ufl", "ufl=2019.1.0", "conda-forge")

function __init__()
	@require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("jplot.jl")
	@require ProgressMeter="92933f4c-e287-5a05-a399-4b506db050ca" begin
	  using ProgressMeter
	end
end
#the below code is an adaptation of aleadev.FEniCS.jl
import Base: size, length, show, *, +, -,/, repr, div, sqrt,split,write
abstract type
  fenicsobject
end #creates placeholder for the fenicsobject type

fenicspycall(object::fenicsobject, func::Union{Symbol,String}, args...) = getproperty(object.pyobject,func)(args...)
export fenicspycall

macro fenicsclass(name::Symbol, base1::Symbol=:fenicsobject)
  impl = Symbol(name, "Impl")
  esc(quote
    abstract type
       $name <: $base1
     end
    struct $impl <: $name
      pyobject::PyObject
    end
    $(name)(pyobject::PyObject) = $impl(pyobject)
  end)
end
export fenicsclass

str(obj::fenicsobject) = fenicspycall(obj, :__str__)
repr(obj::fenicsobject) = fenicspycall(obj, :__repr__)
show(io::IO, obj::fenicsobject) = show(io, str(obj))
Docs.getdoc(obj::fenicsobject) = obj.pyobject.__doc__
export str, repr

include("jmesh.jl") #this file contains the mesh functions
include("jfem.jl") #this file contains the fem functions
include("jmisc.jl") #this file contains various miscallaneous functions to assist with solving etc
include("jsolve.jl") #this file contains the solver functions/routines
include("jinterface.jl")

try
  pyimport_conda("mshr", "mshr=2019.1.0", "conda-forge")
  include("fmshr.jl") #this file contains various geometrical objects using the mshr package
catch ee
 print("mshr has not been included")
end


end #module
