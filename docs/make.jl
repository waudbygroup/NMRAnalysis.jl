using Pkg
Pkg.activate(@__DIR__)  # Activate the docs environment
Pkg.instantiate()       # Install all dependencies

using Documenter, NMRAnalysis

ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

# DocMeta.setdocmeta!(NMRTools, :DocTestSetup, :(using NMRAnalysis); recursive=true)

makedocs(;
         modules=[NMRAnalysis],
         format=Documenter.HTML(),
         pages=["Home" => "index.md",
                "Quick Start" => "quickstart.md",
                "Analyses" => ["Diffusion (1D)" => "analyses/diffusion.md",
                               "TRACT (1D)" => "analyses/tract.md",
                               "R1ρ (1D)" => "analyses/r1rho.md",
                               "2D experiments" => "analyses/2d-fitting.md"],
                "Tutorials" => ["19F R1ρ acquisition & analysis" => "tutorials/r1rho.md"],
                "Advanced" => ["Creating new 2D analyses" => "advanced/creating_2d_analyses.md",
                               "API" => "api.md",
                               "Index" => "indexes.md"]],
         sitename="NMRAnalysis.jl",
         authors="Chris Waudby",
         warnonly=[:missing_docs],)

deploydocs(;
           repo="github.com/waudbygroup/NMRAnalysis.jl.git",
           devbranch="main",
           versions=["stable" => "v^", "v#.#", "dev" => "main"],)
