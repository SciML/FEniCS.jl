using SciMLTesting

run_tests(;
    # Core orchestrates the example sweep and inline FEniCS testsets in a single
    # body (core.jl) rather than as independent per-file groups, so it is passed
    # explicitly instead of being discovered as a folder.
    core = joinpath(@__DIR__, "core.jl"),
    groups = Dict(
        "QA" => (;
            env = joinpath(@__DIR__, "qa"),
            body = joinpath(@__DIR__, "qa", "qa.jl"),
        ),
    ),
    # The original ran only the Core body for the default GROUP="All"; QA ran only
    # when explicitly selected (it exited early before loading FEniCS).
    all = ["Core"],
)
