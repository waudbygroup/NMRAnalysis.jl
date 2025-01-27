"""
    MaybeVector{T} <: AbstractVector{T}

Abstract type representing a vector that may contain either a single element or multiple elements.
All indices return the same value for SingleElementVector.
"""
abstract type MaybeVector{T} <: AbstractVector{T} end

"""
    SingleElementVector{T} <: MaybeVector{T}

Vector type that returns its single element for any valid index.
"""
struct SingleElementVector{T} <: MaybeVector{T}
    x::Vector{T}
    function SingleElementVector{T}(x::Vector{T}) where T
        length(x) == 1 || throw(ArgumentError("SingleElementVector must contain exactly one element"))
        new(x)
    end
end

"""
    StandardVector{T} <: MaybeVector{T}

Standard vector implementation wrapping Vector{T}.
"""
struct StandardVector{T} <: MaybeVector{T}
    x::Vector{T}
end

# Constructors
MaybeVector(x::T) where T = SingleElementVector{T}([x])
MaybeVector(x::Vector{T}) where T = StandardVector{T}(x)

# SingleElementVector interface
Base.size(mv::SingleElementVector) = (1,)
Base.getindex(mv::SingleElementVector, i::Int) = mv.x[1]
Base.setindex!(mv::SingleElementVector, v, i::Int) = (mv.x[1] = v)
Base.length(mv::SingleElementVector) = 1
Base.IndexStyle(::Type{<:SingleElementVector}) = IndexLinear()
Base.similar(mv::SingleElementVector{T}) where T = SingleElementVector{T}(similar(mv.x))
Base.copy(mv::SingleElementVector) = SingleElementVector(copy(mv.x))

# StandardVector interface
Base.size(mv::StandardVector) = size(mv.x)
Base.getindex(mv::StandardVector, i::Int) = mv.x[i]
Base.setindex!(mv::StandardVector, v, i::Int) = (mv.x[i] = v)
Base.length(mv::StandardVector) = length(mv.x)
Base.IndexStyle(::Type{<:StandardVector}) = IndexLinear()
Base.similar(mv::StandardVector{T}) where T = StandardVector{T}(similar(mv.x))
Base.copy(mv::StandardVector) = StandardVector(copy(mv.x))

# Conversions
Base.convert(::Type{StandardVector{T}}, x::SingleElementVector{T}) where T = 
    StandardVector(Vector{T}([x.x[1]]))