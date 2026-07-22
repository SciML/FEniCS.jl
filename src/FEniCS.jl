"""
    FEniCS

Julia wrappers around FEniCS/DOLFIN finite element meshes, function spaces,
forms, solvers, and related helper utilities.
"""
module FEniCS

using PyCall: PyCall, PyObject, pyimport_conda, pycall, @pyimport
using Requires: @require

# determine if we can include mshr
include_mshr = false
try
    pyimport_conda("mshr", "mshr=2019.1.0", "conda-forge")
    global include_mshr = true
catch ee
    print("mshr has not been included")
end

# global PyObject constants that get initialized at runtime.  We
# initialize them here (rather than via "global foo = ..." in __init__)
# so that their type is known at compile-time.
const fenics = PyCall.PyNULL()
const ufl = PyCall.PyNULL()
const mshr = PyCall.PyNULL()

function __init__()
    @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" include("jplot.jl")
    copy!(fenics, pyimport_conda("fenics", "fenics=2019.1.0", "conda-forge"))
    copy!(ufl, pyimport_conda("ufl", "ufl=2019.1.0", "conda-forge"))
    include_mshr && copy!(mshr, pyimport_conda("mshr", "mshr=2019.1.0", "conda-forge"))

    # initialize constants at loadtime
    global dx = Measure(fenics.dx)
    global ds = Measure(fenics.ds)
    global dS = Measure(fenics.dS)
    global dP = Measure(fenics.dP)
    global tetrahedron = fenics.tetrahedron
    global hexahedron = fenics.hexahedron #matplotlib cannot handle hexahedron elements
    global triangle = fenics.triangle
    global quadrilateral = fenics.quadrilateral
    return copy!(CellType.pyobject, fenics.CellType)
end
#the below code is an adaptation of aleadev.FEniCS.jl
import Base: size, length, show, *, +, -, /, ^, sin, cos, tan, asin, acos, atan, exp, log,
    repr, div, sqrt, split, write

import SpecialFunctions: besseli, besselj, besselk, bessely
export besseli, besselj, besselk, bessely

import LinearAlgebra: norm
export norm

abstract type fenicsobject end

struct PyNamespace <: fenicsobject
    pyobject::PyObject
end

function Base.getproperty(obj::PyNamespace, name::Symbol)
    name === :pyobject && return getfield(obj, :pyobject)
    return getproperty(getfield(obj, :pyobject), name)
end

"""
    CellType

FEniCS/DOLFIN cell-type namespace used to select mesh cell shapes.
"""
const CellType = PyNamespace(PyCall.PyNULL())

function fenicspycall(object::fenicsobject, func::Union{Symbol, String}, args...)
    return getproperty(object.pyobject, func)(args...)
end
export fenicspycall

"""
    @fenicsclass name [base = FEniCS.fenicsobject]

Define a Julia wrapper hierarchy for a FEniCS Python object.

# Arguments
- `name`: Name of the abstract wrapper type to define.
- `base`: Optional abstract supertype. It defaults to FEniCS's wrapper base type.

# Extension rules
The macro defines `name <: base`, a concrete `nameImpl <: name` with one
`PyCall.PyObject` field named `pyobject`, and a constructor from that field.
Use it from a module that imports `FEniCS`; importing `PyCall` is not required
to expand the macro. Extend wrapper behavior on `name`, not on `nameImpl`, so
all implementations of the abstract wrapper share the extension.

# Example
```julia
module MyFEniCSExtension
using FEniCS
FEniCS.@fenicsclass MyObject
end
```
"""
macro fenicsclass(name::Symbol, base1 = nothing)
    impl = Symbol(name, "Impl")
    base = base1 === nothing ? GlobalRef(@__MODULE__, :fenicsobject) : esc(base1)
    pyobject = GlobalRef(PyCall, :PyObject)
    return quote
        abstract type $(esc(name)) <: $base end
        struct $(esc(impl)) <: $(esc(name))
            pyobject::$pyobject
        end
        $(esc(name))(pyobject::$pyobject) = $(esc(impl))(pyobject)
    end
end
export @fenicsclass

@enum LOGLEVEL begin
    CRITICAL = 50
    ERROR = 40
    WARNING = 30
    INFO = 20
    PROGRESS = 16
    TRACE = 13
    DBG = 10
end
set_log_level(lvl::LOGLEVEL) = set_log_level(Int(lvl))
set_log_level(lvl::Int) = fenics.set_log_level(lvl)

str(obj::fenicsobject) = fenicspycall(obj, :__str__)
repr(obj::fenicsobject) = fenicspycall(obj, :__repr__)
show(io::IO, obj::fenicsobject) = show(io, repr(obj))
Docs.getdoc(::PyNamespace) = "FEniCS/DOLFIN cell-type namespace used to select mesh cell shapes."
function Docs.getdoc(obj::fenicsobject)
    doc = convert(PyCall.PyAny, obj.pyobject.__doc__)
    return doc isa AbstractString ? doc : "FEniCS/DOLFIN wrapped object."
end
export str, repr

include("jmesh.jl") #this file contains the mesh functions
include("jfem.jl") #this file contains the fem functions
include("jmisc.jl") #this file contains various miscallaneous functions to assist with solving etc
include("jsolve.jl") #this file contains the solver functions/routines
include("jinterface.jl")
include_mshr && include("fmshr.jl") #this file contains various geometrical objects using the mshr package
include("public_api_docs.jl")

end #module
