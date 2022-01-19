using GrayCoding
using Documenter

DocMeta.setdocmeta!(GrayCoding, :DocTestSetup, :(using GrayCoding); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/nivupai/GrayCoding.jl",
    devbranch="main",
)
