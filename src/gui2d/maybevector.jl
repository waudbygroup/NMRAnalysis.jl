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
mutable struct SingleElementVector{T} <: MaybeVector{T}
    x::T
end

"""
    StandardVector{T} <: MaybeVector{T}

Standard vector implementation wrapping Vector{T}.
"""
struct StandardVector{T} <: MaybeVector{T}
    x::Vector{T}
end

# Constructors
MaybeVector(x::T) where T = SingleElementVector{T}(x)
MaybeVector(x::Vector{T}) where T = StandardVector{T}(x)

# SingleElementVector interface
Base.size(mv::SingleElementVector) = (1,)
Base.getindex(mv::SingleElementVector, i::Int) = mv.x
Base.setindex!(mv::SingleElementVector, v, i::Int) = (mv.x = v)
Base.length(mv::SingleElementVector) = 1
Base.IndexStyle(::Type{<:SingleElementVector}) = IndexLinear()
Base.similar(mv::SingleElementVector{T}) where T = SingleElementVector{T}(mv.x)
Base.copy(mv::SingleElementVector) = SingleElementVector(mv.x)

# StandardVector interface
Base.size(mv::StandardVector) = size(mv.x)
Base.getindex(mv::StandardVector, i::Int) = mv.x[i]
Base.setindex!(mv::StandardVector, v, i::Int) = (mv.x[i] = v)
Base.length(mv::StandardVector) = length(mv.x)
Base.IndexStyle(::Type{<:StandardVector}) = IndexLinear()
Base.similar(mv::StandardVector{T}) where T = StandardVector{T}(similar(mv.x))
Base.copy(mv::StandardVector) = StandardVector(copy(mv.x))

# Iteration protocol
Base.iterate(mv::SingleElementVector) = (mv.x, nothing)
Base.iterate(mv::SingleElementVector, state) = nothing
Base.iterate(mv::StandardVector) = iterate(mv.x)
Base.iterate(mv::StandardVector, state) = iterate(mv.x, state)

# Broadcasting
Base.broadcastable(mv::SingleElementVector) = Ref(mv.x)
Base.broadcastable(mv::StandardVector) = mv.x

# Conversions
Base.convert(::Type{StandardVector{T}}, x::SingleElementVector{T}) where T = 
    StandardVector([x.x])
Base.convert(::Type{SingleElementVector{T}}, x::StandardVector{T}) where T = 
    length(x.x) == 1 ? SingleElementVector(first(x.x)) : 
    throw(ArgumentError("Cannot convert StandardVector with length != 1 to SingleElementVector"))
Base.convert(::Type{Vector{T}}, x::SingleElementVector{T}) where T = [x.x]
Base.convert(::Type{Vector{T}}, x::StandardVector{T}) where T = x.x
Base.convert(::Type{MaybeVector{T}}, x::Vector{T}) where T = 
    length(x) == 1 ? SingleElementVector(first(x)) : StandardVector(x)


# # Unit Tests
# using Test

# @testset "MaybeVector Tests" begin
#     @testset "Construction" begin
#         s = MaybeVector(1)
#         @test s isa SingleElementVector{Int}
#         @test s.x == 1
        
#         v = MaybeVector([1, 2, 3])
#         @test v isa StandardVector{Int}
#         @test v.x == [1, 2, 3]
#     end

#     @testset "Indexing" begin
#         s = MaybeVector(1)
#         @test s[1] == 1
#         @test s[2] == 1  # All valid indices return the same value
        
#         v = MaybeVector([1, 2, 3])
#         @test v[1] == 1
#         @test v[2] == 2
#         @test v[3] == 3
#         @test_throws BoundsError v[4]
#     end

#     @testset "Length and Size" begin
#         s = MaybeVector(1)
#         @test length(s) == 1
#         @test size(s) == (1,)
        
#         v = MaybeVector([1, 2, 3])
#         @test length(v) == 3
#         @test size(v) == (3,)
#     end

#     @testset "Iteration" begin
#         s = MaybeVector(1)
#         collected_s = collect(s)
#         @test collected_s == [1]
        
#         v = MaybeVector([1, 2, 3])
#         collected_v = collect(v)
#         @test collected_v == [1, 2, 3]
#     end

#     @testset "Broadcasting" begin
#         s = MaybeVector(2)
#         @test s .* 2 == 4
        
#         v = MaybeVector([1, 2, 3])
#         @test v .* 2 == [2, 4, 6]
#     end

#     @testset "Conversion" begin
#         # Single -> Standard
#         s = MaybeVector(1)
#         converted_to_standard = convert(StandardVector{Int}, s)
#         @test converted_to_standard isa StandardVector{Int}
#         @test converted_to_standard.x == [1]
        
#         # Standard -> Single (valid case)
#         v = MaybeVector([1])
#         converted_to_single = convert(SingleElementVector{Int}, v)
#         @test converted_to_single isa SingleElementVector{Int}
#         @test converted_to_single.x == 1
        
#         # Standard -> Single (invalid case)
#         v_long = MaybeVector([1, 2])
#         @test_throws ArgumentError convert(SingleElementVector{Int}, v_long)
        
#         # To Vector
#         @test convert(Vector{Int}, s) == [1]
#         @test convert(Vector{Int}, v) == [1]
        
#         # From Vector
#         @test convert(MaybeVector{Int}, [1]) isa SingleElementVector{Int}
#         @test convert(MaybeVector{Int}, [1, 2]) isa StandardVector{Int}
#     end

#     @testset "Mutation" begin
#         s = MaybeVector(1)
#         s[1] = 2
#         @test s[1] == 2
#         @test s[2] == 2  # All indices show the new value
        
#         v = MaybeVector([1, 2, 3])
#         v[2] = 4
#         @test v[2] == 4
#         @test v.x == [1, 4, 3]
#     end
# end