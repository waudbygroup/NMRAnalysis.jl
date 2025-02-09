using Documenter#, NMRAnalysis
ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

# DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRAnalysis); recursive=true)

makedocs(;
    modules=[NMRAnalysis],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Reference guide" => [
            "Creating new experiments" => "experiments.md",
        ],
        "API" => "api.md",
        "Index" => "indexes.md",
    ],
    sitename="NMRAnalysis.jl",
    authors="Chris Waudby",
    assets=String[],
    warnonly = [:missing_docs],
)

deploydocs(;
    repo="github.com/waudbygroup/NMRAnalysis.jl.git",
    versions = ["stable" => "v^", "v#.#", "dev" => "main"],
)