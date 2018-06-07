__precompile__(false)
module FEniCS
using PyCall
@pyimport fenics

#the below code is an adaptation of aleadev.FEniCS.jl
import Base: size, length, show, *, +, -,/, repr, dot, cross, div,Function, sqrt
abstract type
  fenicsobject
end #creates placeholder for the fenicsobject type

fenicspycall(object::fenicsobject, func::Union{Symbol,String}, args...) = object.pyobject[func](args...)
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
export str, repr

include("jmesh.jl") #this file contains the mesh functions
include("jfem.jl") #this file contains the fem functions
include("jmisc.jl") #this file contains various miscallaneous functions to assist with solving etc
include("jsolve.jl") #this file contains the solver functions/routines
include("jinterface.jl")
try
  pyimport("mshr")
  include("fmshr.jl") #this file contains various geometrical objects using the mshr package
catch ee
 print("mshr has not been included")
end
try
  Pkg.installed("PyPlot")
  include("jplot.jl")
  print("plotting loaded")
catch ee
  print("PyPlot is not installed. plotting is not available")
end

end #module
