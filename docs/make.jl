using Polylogarithms
using Documenter

DocMeta.setdocmeta!(Polylogarithms, :DocTestSetup, :(using Polylogarithms); recursive=true)

makedocs(;
    modules=[Polylogarithms],
    authors="Matthew Roughan <matthew.roughan@adelaide.edu.au>",
    sitename="Polylogarithms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://matthew.roughan@adelaide.edu.au.github.io/Polylogarithms.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mroughan/Polylogarithms.jl",
)
