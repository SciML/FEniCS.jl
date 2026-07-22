using SciMLTesting, FEniCS, Test
using JET

const _DEPENDENCY_REEXPORTS = (
    :besseli,
    :besselj,
    :besselk,
    :bessely,
    :div,
    :norm,
    :repr,
    :size,
    :split,
    :sqrt,
    :write,
)

run_qa(
    FEniCS;
    reexports_allow = _DEPENDENCY_REEXPORTS,
    api_docs_kwargs = (; rendered_ignore = _DEPENDENCY_REEXPORTS),
    ei_kwargs = (;
        # `getdoc` is not public in `Base.Docs`, but FEniCS extends
        # `Base.Docs.getdoc(::fenicsobject)` (src/FEniCS.jl) to surface the wrapped
        # PyObject's docstring. Extending a Base-internal method legitimately
        # requires naming it; nothing to make public here.
        all_qualified_accesses_are_public = (; ignore = (:getdoc,)),
    ),
)
