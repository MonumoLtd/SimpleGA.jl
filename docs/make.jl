using GeometricAlgebra
using Documenter

DocMeta.setdocmeta!(
    GeometricAlgebra, :DocTestSetup, :(using GeometricAlgebra); recursive=true
)

makedocs(;
    modules=[GeometricAlgebra],
    authors="Monumo Ltd.",
    repo="https://github.com/MonumoLtd/GeometricAlgebra.jl/blob/{commit}{path}#{line}",
    sitename="GeometricAlgebra.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://monumoltd.github.io/GeometricAlgebra.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=["Home" => "index.md"],
    checkdocs=:exports,
    strict=true,
)

deploydocs(; repo="github.com/MonumoLtd/GeometricAlgebra.jl", devbranch="main")
