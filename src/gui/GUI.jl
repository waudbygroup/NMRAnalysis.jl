module GUI

using GLMakie
using NMRTools

export view2d
export fit2d

include("util.jl")
include("types.jl")
include("peaks.jl")
include("clustering.jl")
include("specdata.jl")
include("view2d.jl")
include("fit2d.jl")

end