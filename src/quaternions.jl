"""
    Quaternions

Quaternion code. This is called by some of the later GA implementations.

The core mirrors much of the GA code structure.
For completeness we have defined a division operation for quaternions as they are a division
algebra.
"""
module Quaternions

using SimpleGA
using LinearAlgebra

export real_part, imag_part, Quaternion

struct Quaternion{T<:Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

Quaternion{T}(w::Real) where {T<:Real} = Quaternion(convert(T, w))
Quaternion(w::Real) = Quaternion(w, zero(w), zero(w), zero(w))
Quaternion(w::Real, x::Real, y::Real, z::Real) = Quaternion(promote(w, x, y, z)...)

function Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real,S<:Real}
    return Quaternion{promote_type(T, S)}
end
function Base.promote_rule(
    ::Type{Quaternion{T}}, ::Type{Quaternion{S}}
) where {T<:Real,S<:Real}
    return Quaternion{promote_type(T, S)}
end

function Base.convert(::Type{Quaternion{T}}, a::Quaternion) where {T<:Real}
    return Quaternion{T}(convert(T, a.w), convert(T, a.x), convert(T, a.y), convert(T, a.z))
end

Base.zero(a::Quaternion) = Quaternion(zero(a.w))
Base.one(a::Quaternion) = Quaternion(one(a.w))

Base.:(-)(a::Quaternion) = Quaternion(-a.w, -a.x, -a.y, -a.z)
function Base.:(+)(a::Quaternion, b::Quaternion)
    return Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end
function Base.:(-)(a::Quaternion, b::Quaternion)
    return Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end
Base.:(+)(x::Number, a::Quaternion) = Quaternion(x + a.w, a.x, a.y, a.z)
Base.:(-)(x::Number, a::Quaternion) = Quaternion(x - a.w, -a.x, -a.y, -a.z)
Base.:(+)(a::Quaternion, x::Number) = Quaternion(a.w + x, a.x, a.y, a.z)
Base.:(-)(a::Quaternion, x::Number) = Quaternion(a.w - x, a.x, a.y, a.z)

#Multiplication
Base.:(*)(x::Number, a::Quaternion) = Quaternion(x * a.w, x * a.x, x * a.y, x * a.z)
Base.:(*)(a::Quaternion, x::Number) = x * a

function Base.:(*)(a::Quaternion, b::Quaternion)
    return Quaternion(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

Base.:(/)(a::Quaternion, x::Number) = (1 / x) * a
Base.:(/)(a::Quaternion, b::Quaternion) = a * conj(b) / dot(b, conj(b))
LinearAlgebra.norm(a::Quaternion) = sqrt(a.w^2 + a.x^2 + a.y^2 + a.z^2)

#Reverse
Base.conj(a::Quaternion) = Quaternion(a.w, -a.x, -a.y, -a.z)

#Projection operations
function LinearAlgebra.dot(a::Quaternion, b::Quaternion)
    return a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
end
Base.real(a::Quaternion) = a.w
real_part(a::Quaternion) = Quaternion(a.w, zero(a.w), zero(a.w), zero(a.w))
imag_part(a::Quaternion) = Quaternion(zero(a.w), a.x, a.y, a.z)

#Exponentiation
function SimpleGA.bivector_exp(a::Quaternion)
    a = imag_part(a)
    nrm = norm(a)
    return if iszero(nrm)
        Quaternion(one(nrm), 0, 0, 0)
    else
        cos(nrm) + sin(nrm) * a / nrm
    end
end

function Base.exp(a::Quaternion)
    R = SimpleGA.bivector_exp(a)
    return iszero(a.w) ? R : exp(a.w) * R
end

#Additional Functions
function Base.isapprox(a::Quaternion, b::Quaternion)
    return isapprox(a.w, b.w) &&
           isapprox(a.x, b.x) &&
           isapprox(a.y, b.y) &&
           isapprox(a.z, b.z)
end

function Base.isapprox(a::Quaternion{T}, b::Quaternion{T}; kwargs...) where {T}
    return (
        isapprox(a.w, b.w; kwargs...) &&
        isapprox(a.x, b.x; kwargs...) &&
        isapprox(a.y, b.y; kwargs...) &&
        isapprox(a.z, b.z; kwargs...)
    )
end

function Base.show(io::IO, ::MIME"text/plain", a::Quaternion)
    return print(
        io,
        "",
        string(a.w) *
        " + " *
        string(a.x) *
        "i + " *
        string(a.y) *
        "j + " *
        string(a.z) *
        "k",
    )
end

function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Quaternion})
    n = length(mvs)
    println(io, n, "-element Vector{Quaternion}")
    for i in eachindex(mvs)
        println(io, " ", mvs[i])
    end
end

end #Module
