using GADraft
using Documenter

DocMeta.setdocmeta!(GADraft, :DocTestSetup, :(using GADraft); recursive=true)

makedocs(;
    modules=[GADraft],
    authors="Tom Gillam <tpgillam@googlemail.com>",
    repo="https://github.com/tpgillam/GADraft.jl/blob/{commit}{path}#{line}",
    sitename="GADraft.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tpgillam.github.io/GADraft.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    strict=true,
)

deploydocs(;
    repo="github.com/tpgillam/GADraft.jl",
    devbranch="main",
)
