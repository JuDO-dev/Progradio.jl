using Progradio
using Documenter

DocMeta.setdocmeta!(Progradio, :DocTestSetup, :(using Progradio); recursive=true)

makedocs(;
    modules=[Progradio],
    authors="astroEduardo <72969764+astroEduardo@users.noreply.github.com> and contributors",
    repo="https://github.com/JuDO-dev/Progradio.jl/blob/{commit}{path}#{line}",
    sitename="Progradio.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuDO-dev.github.io/Progradio.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuDO-dev/Progradio.jl",
    devbranch="dev",
)
