using Documenter
ENV["GKSwstype"] = "100" # https://github.com/jheinen/GR.jl/issues/278

makedocs(;
         modules=[NMRAnalysis],
         format=Documenter.HTML(),
         pages=["Home" => "index.md",
                "R1ρ" => "r1rho.md",
                "Tutorials" => ["R1ρ analysis" => "tutorials/R1rho.md"]],
         sitename="NMRAnalysis.jl",
         authors="Chris Waudby",
         assets=String[],
         warnonly=[:missing_docs],)

deploydocs(;
           repo="github.com/waudbygroup/NMRAnalysis.jl.git",
           versions=["stable" => "v^", "v#.#", "dev" => "main"],)
