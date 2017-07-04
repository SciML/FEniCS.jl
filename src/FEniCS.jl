
module FEniCS
using PyCall
@pyimport fenics
#the below code is an adaptation of aleadev.FEniCS.jl

import Base: size, length, show, *, +, -, Function,repr
abstract fenicsobject #creates placeholder for the fenicsobject type
fenicspycall(object::fenicsobject, func::Union{Symbol,String}, args...) = object.pyobject[func](args...)

macro fenicsclass(name::Symbol, base1::Symbol=:fenicsobject, base2::Symbol=:Any, base3::Symbol=:Any)
  impl = Symbol(name, "Impl")
  esc(quote
    abstract $name <: $base1, $base2, $base3
    immutable $impl <: $name
      pyobject::PyObject
    end
    $(name)(pyobject::PyObject) = $impl(pyobject)
  end)
end

str(obj::fenicsobject) = fenicspycall(obj, :__str__)
repr(obj::fenicsobject) = fenicspycall(obj, :__repr__)
show(io::IO, obj::fenicsobject) = show(io, str(obj))
export str, repr
#the code from aleadev.fenics.jl has ended. The code underneath should be moved to the mesh.jl file before being uploaded
include("jmesh.jl")
end# module
