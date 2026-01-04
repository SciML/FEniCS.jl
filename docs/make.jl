using Documenter

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

makedocs(
    sitename = "FEniCS.jl",
    authors = "Chris Rackauckas",
    modules = Module[],
    clean = true, doctest = false, linkcheck = true,
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
