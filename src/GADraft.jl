module GADraft

# Core types and functions.
export Parity, Even, Odd, MultiVector
export parity, iseven, isodd
export project, scalar

# Concrete Geometric algebra modules
export GA20, GA30

"""
    Parity

A concrete [`MultiVector`](@ref) will represent either even grade terms, or odd grade terms,
but never a mixture.

The `Parity` abstract type is a supertype of [`Even`](@ref) and [`Odd`](@ref).
"""
abstract type Parity end

struct Even <: Parity end
struct Odd <: Parity end

"""
    MultiVector{P<:Parity,T<:Real}

TODO
TODO
TODO
"""
abstract type MultiVector{P<:Parity,T<:Real} end

Base.iseven(::MultiVector{Even}) = true
Base.iseven(::MultiVector{Odd}) = false
Base.isodd(mv::MultiVector) = !iseven(mv)

# TODO Check if we actually use these anywhere, delete if not.
# Allow inversion of type parameters for convenience.
Base.:(!)(::Type{Even}) = Odd
Base.:(!)(::Type{Odd}) = Even

# For internal use only (not exported).
_scalar_T(::MultiVector{P,T}) where {P,T} = T
_complex_T(::Complex{T}) where {T<:Real} = T

# Used in implementations of Base.show.
function string_parts(x::Real, label::Union{Nothing,AbstractString}=nothing)
    return if iszero(x)
        String[]
    elseif isnothing(label)
        ["$x"]
    else
        ["$x $label"]
    end
end

"""
    project(mv::MultiVector, grade::Integer) -> MultiVector

Project out the part of the multivector `mv` of the given `grade`.
"""
function project end

"""
    scalar(mv::MultiVector) -> Real
    scalar(a::MultiVector, b::MultiVector) -> Real

In the first form, return the scalar part of `mv`.

In the second form, return the scalar part of `a * b`. Whilst equivalent to `scalar(a * b)`,
`scalar(a, b)` avoids unnecessary computation, so will be more performant.
"""
function scalar end

# Implementations that we know will give zero regardless of dimension.
scalar(a::MultiVector{Odd}) = zero(_scalar_T(a))
function scalar(a::MultiVector{Odd}, b::MultiVector{Even})
    return zero(promote_type(_scalar_T(a), _scalar_T(b)))
end
function scalar(a::MultiVector{Even}, b::MultiVector{Odd})
    return zero(promote_type(_scalar_T(a), _scalar_T(b)))
end

# Lightweight quaternion implementation, which is used in certain algebras.
include("quaternion.jl")

# Implementations of algebras.
include("ga20.jl")
include("ga30.jl")

end
