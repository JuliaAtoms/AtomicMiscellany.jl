using AtomicMiscellany
using Documenter

DocMeta.setdocmeta!(AtomicMiscellany, :DocTestSetup, quote
    using AtomicMiscellany, AtomicMiscellany.NuclearModels
end; recursive=true)

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
        "Implementation notes" => [
            "implementation/nuclearmodels.md",
        ]
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/AtomicMiscellany.jl",
)
