using Polylogarithms
using Documenter

makedocs(;
    modules=[Polylogarithms],
    authors="Matthew Roughan <matthew.roughan@adelaide.edu.au> and contributors",
    repo="https://github.com/matthew.roughan@adelaide.edu.au/Polylogarithms.jl/blob/{commit}{path}#L{line}",
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
    repo="github.com/matthew.roughan@adelaide.edu.au/Polylogarithms.jl",
)
