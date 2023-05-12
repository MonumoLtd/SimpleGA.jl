using SimpleGA
using Documenter

DocMeta.setdocmeta!(SimpleGA, :DocTestSetup, :(using SimpleGA); recursive=true)

makedocs(;
    modules=[SimpleGA],
    authors="Monumo Ltd.",
    repo="https://github.com/MonumoLtd/SimpleGA.jl/blob/{commit}{path}#{line}",
    sitename="SimpleGA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://monumoltd.github.io/SimpleGA.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md", "Reference" => "reference.md"],
    checkdocs=:exports,
    strict=true,
)

deploydocs(; repo="github.com/MonumoLtd/SimpleGA.jl", devbranch="main")
