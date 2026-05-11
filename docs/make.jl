using Documenter

makedocs(
    sitename = "FastEarth.jl",
    authors  = "Alexander Robinson and contributors",
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        mathengine = Documenter.MathJax3(),
        canonical  = "https://fesmc.github.io/FastEarth.jl/",
        edit_link  = "main",
    ),
    pages = [
        "Home"   => "index.md",
        "Design" => "design.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo      = "github.com/fesmc/FastEarth.jl.git",
    devbranch = "main",
    push_preview = true,
)
