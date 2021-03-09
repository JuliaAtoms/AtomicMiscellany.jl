using AtomicMiscellany
using Documenter

DocMeta.setdocmeta!(AtomicMiscellany, :DocTestSetup, :(using AtomicMiscellany); recursive=true)

makedocs(;
    modules=[AtomicMiscellany],
    authors="JuliaAtoms contributors",
    repo="https://github.com/JuliaAtoms/AtomicMiscellany.jl/blob/{commit}{path}#{line}",
    sitename="AtomicMiscellany.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaAtoms.github.io/AtomicMiscellany.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/AtomicMiscellany.jl",
)
