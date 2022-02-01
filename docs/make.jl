using GrayCoding
using Documenter

DocMeta.setdocmeta!(GrayCoding, :DocTestSetup, :(using GrayCoding); recursive=true)
push!(LOAD_PATH,"../src/")
makedocs(;
    modules=[GrayCoding],
    authors="Nivedita Rethnakar et al.",
    repo="https://github.com/nivupai/GrayCoding.jl/blob/{commit}{path}#{line}",
    sitename="GrayCoding.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nivupai.github.io/GrayCoding.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Algebra of Gray Codes" => "algebra.md",
        "Applications" => [
            "List of Applications" => "applications.md",
            "Quantum Algorithms and Circuits" => "quantum.md"
        ],
        "Tutorials" => "tutorials.md",
    ],
)

deploydocs(;
    repo="github.com/nivupai/GrayCoding.jl",
    devbranch="main",
)
