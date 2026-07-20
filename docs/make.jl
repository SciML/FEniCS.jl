using Documenter

try
    using FEniCS
catch err
    @warn "FEniCS failed to initialize; rendering static API docstrings." exception = err

    @eval module FEniCS
    """
        @fenicsclass name [base]

    Define a Julia wrapper type hierarchy for a FEniCS Python object.
    """
    macro fenicsclass(args...)
        return nothing
    end

    include("../src/public_api_docs.jl")

    const _DOCS_BUILD_API_DOCS = Pair{Symbol, String}[
        :Box => "Create a three-dimensional mshr box geometry from two opposite corners.",
        :CSGDifference => "Construct the difference of two mshr constructive-solid-geometry objects.",
        :CSGIntersection => "Construct the intersection of two mshr constructive-solid-geometry objects.",
        :CSGRotation => "Construct a rotated mshr constructive-solid-geometry object.",
        :CSGScaling => "Construct a scaled mshr constructive-solid-geometry object.",
        :CSGTranslation => "Construct a translated mshr constructive-solid-geometry object.",
        :CSGUnion => "Construct the union of two mshr constructive-solid-geometry objects.",
        :CellType => "FEniCS/DOLFIN cell-type namespace used to select mesh cell shapes.",
        :Circle => "Create a two-dimensional mshr circle geometry.",
        :Cone => "Create a three-dimensional mshr cone geometry.",
        :Ellipse => "Create a two-dimensional mshr ellipse geometry.",
        :Extrude2D => "Extrude a two-dimensional mshr geometry into three dimensions.",
        :Rectangle => "Create a two-dimensional mshr rectangle geometry.",
        :Sphere => "Create a three-dimensional mshr sphere geometry.",
        :feMesh => "Mesh data container used by the Julia finite element interface.",
        :generate_mesh => "Generate a FEniCS mesh from an mshr geometry object.",
        :plot => "Plot a FEniCS mesh, expression, function, or finite-element solution.",
        :set_subdomain => "Set a subdomain on an mshr geometry object.",
        :surf_plot => "Create a surface plot for a scalar FEniCS expression or function.",
    ]

    end

    for (name, doc) in FEniCS._DOCS_BUILD_API_DOCS
        Core.eval(FEniCS, :(@doc $doc $name))
    end

    for (name, _) in Iterators.flatten((FEniCS._PUBLIC_API_DOCS, FEniCS._DOCS_BUILD_API_DOCS))
        isdefined(FEniCS, name) ||
            Core.eval(FEniCS, Expr(:const, Expr(:(=), name, nothing)))
    end
end

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

makedocs(
    sitename = "FEniCS.jl",
    authors = "Chris Rackauckas",
    modules = Module[],
    clean = true, doctest = false, linkcheck = true,
    linkcheck_ignore = [
        r"https://github\.com/.*",
        r"https://gist\.github\.com/.*",
    ],
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"],
        canonical = "https://docs.sciml.ai/FEniCS/stable/"
    ),
    pages = [
        "FEniCS.jl: Finite Element PDE Solving in Julia" => "index.md",
    ]
)

deploydocs(;
    repo = "github.com/SciML/FEniCS.jl"
)
