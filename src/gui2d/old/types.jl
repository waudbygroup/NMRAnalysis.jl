abstract type Model end




abstract type MaybeVector{T} <: AbstractVector{T} end

struct SingleElementVector{T} <: MaybeVector{T}
    x::Vector{T}
end

struct StandardVector{T} <: MaybeVector{T}
    x::Vector{T}
end

# Constructor for MaybeVector
MaybeVector(x::T) where T = SingleElementVector{T}([x])
MaybeVector(x::Vector{T}) where T = StandardVector{T}(x)

# Implement AbstractVector interface for SingleElementVector
Base.size(mv::SingleElementVector) = (1,)
Base.getindex(mv::SingleElementVector, i::Int) = mv.x[1]
Base.setindex!(mv::SingleElementVector, v, i::Int) = (mv.x[1] = v)
Base.length(mv::SingleElementVector) = 1
Base.IndexStyle(::Type{SingleElementVector}) = IndexLinear()

# Implement AbstractVector interface for StandardVector
Base.size(mv::StandardVector) = size(mv.x)
Base.getindex(mv::StandardVector, i::Int) = mv.x[i]
Base.setindex!(mv::StandardVector, v, i::Int) = (mv.x[i] = v)
Base.length(mv::StandardVector) = length(mv.x)
Base.IndexStyle(::Type{StandardVector}) = IndexLinear()



mutable struct PeakPseudo2D
    position::Point2f
    amp::Vector{Float64}
    lwX::Float64
    lwY::Float64

    position_err::Point2f
    amp_err::Vector{Float64}
    lwX_err::Float64
    lwY_err::Float64

    A::Float64
    k::Float64
    A_err::Float64
    k_err::Float64

    label::String
    touched::Bool
end
PeakPseudo2D(position::Point2f, amplitude::Vector{Float64}, lwX=30, lwY=30, label="") = PeakPseudo2D(position, amplitude, lwX, lwY, Point2f(0.0, 0.0), 0amplitude, 0.0, 0.0, 0., 0., 0., 0., label, true)
PeakPseudo2D(position::Point2f, label="") = PeakPseudo2D(position, Vector{Float64}(), 30., 30., Point2f(0.0, 0.0), Vector{Float64}(), 0.0, 0.0, 0., 0., 0., 0., label, true)