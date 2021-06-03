using AtomicMiscellany
using Documenter, Literate
using Plots, Formatting

# Stop GR from trying to connect to X when producing plots.
# This will make the documentation build about x10 faster and gets rid of all the
# GR error messages from the build logs.
# https://discourse.julialang.org/t/generation-of-documentation-fails-qt-qpa-xcb-could-not-connect-to-display/60988/2
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(AtomicMiscellany, :DocTestSetup, quote
    using AtomicMiscellany, AtomicMiscellany.NuclearModels
end; recursive=true)

Literate.markdown(
    joinpath(@__DIR__, "src", "implementation", "hydrogenic.jl"),
    joinpath(@__DIR__, "src", "implementation"),
    documenter = true,
)

makedocs(;
    modules=[AtomicMiscellany],
    authors="JuliaAtoms contributors",
    repo="https://github.com/JuliaAtoms/AtomicMiscellany.jl/blob/{commit}{path}#{line}",
    sitename="AtomicMiscellany.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaAtoms.github.io/AtomicMiscellany.jl",
        assets=String[],
        collapselevel = 1,
    ),
    pages=[
        "Home" => "index.md",
        "nuclearmodels.md",
        "hydrogenic.md",
        "Implementation notes" => [
            "Nuclear models" => "implementation/nuclearmodels.md",
            "Hydrogenic energies" => "implementation/hydrogenic.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/AtomicMiscellany.jl",
)
