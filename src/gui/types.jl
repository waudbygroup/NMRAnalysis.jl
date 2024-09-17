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



