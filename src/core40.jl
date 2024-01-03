"""
Core code for the implementation of GA(4,0).

Work using self-dual and anti-self-dual decomposition, so just use a pair of quaternions.
Useful algebra for projective geometry.
"""

using ..Quaternions

struct Even{T<:Real} <: Number
    qp::Quaternion{T}
    qm::Quaternion{T}
end

struct Odd{T<:Real} <: Number
    qp::Quaternion{T}
    qm::Quaternion{T}
end

function Base.convert(::Type{Even{T}}, a::Even) where {T<:Real}
    return Even{T}(convert(Quaternion{T}, a.qp), convert(Quaternion{T}, a.qm))
end
function Base.convert(::Type{Odd{T}}, a::Odd) where {T<:Real}
    return Odd{T}(convert(Quaternion{T}, a.qp), convert(Quaternion{T}, a.qm))
end
function Base.promote_rule(::Type{Even{S}}, ::Type{Even{T}}) where {S<:Real,T<:Real}
    return Even{promote_type(S, T)}
end
function Base.promote_rule(::Type{Odd{S}}, ::Type{Odd{T}}) where {S<:Real,T<:Real}
    return Odd{promote_type(S, T)}
end

Base.zero(a::Even) = Even(zero(a.qp), zero(a.qm))
Base.zero(a::Odd) = Odd(zero(a.qp), zero(a.qm))
Base.one(a::Even) = Even(one(a.qp), one(a.qm))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.qp, -a.qm)
Base.:(-)(a::Odd) = Odd(-a.qp, -a.qm)
Base.:(+)(a::Even, b::Even) = Even(a.qp + b.qp, a.qm + b.qm)
Base.:(+)(a::Odd, b::Odd) = Odd(a.qp + b.qp, a.qm + b.qm)
Base.:(-)(a::Even, b::Even) = Even(a.qp - b.qp, a.qm - b.qm)
Base.:(-)(a::Odd, b::Odd) = Odd(a.qp - b.qp, a.qm - b.qm)

#Scalar addition / subtraction. Other cases are in GAcommon
Base.:(+)(num::Number, a::Even) = Even(a.qp + num, a.qm + num)
Base.:(-)(num::Number, a::Even) = Even(-a.qp + num, -a.qm + num)

#Multiplication
Base.:(*)(num::Number, a::Even) = Even(num * a.qp, num * a.qm)
Base.:(*)(num::Number, a::Odd) = Odd(num * a.qp, num * a.qm)
Base.:(*)(a::Even, b::Even) = Even(a.qp * b.qp, a.qm * b.qm)
Base.:(*)(a::Even, b::Odd) = Odd(a.qp * b.qp, a.qm * b.qm)
Base.:(*)(a::Odd, b::Even) = Odd(a.qp * b.qm, a.qm * b.qp)
Base.:(*)(a::Odd, b::Odd) = Even(a.qp * b.qm, a.qm * b.qp)

#Reverse
LinearAlgebra.adjoint(a::Even) = Even(conj(a.qp), conj(a.qm))
LinearAlgebra.adjoint(a::Odd) = Odd(conj(a.qm), conj(a.qp))

#Grade and projection
function SimpleGA.project(a::Even, n::Integer)
    return if (n == 0)
        tmp = a.qp.w + a.qm.w
        scl = convert(typeof(tmp), tmp / 2)
        scl * one(a)
    elseif (n == 2)
        convert(typeof(a), (a - a') / 2)
    elseif (n == 4)
        tmp = a.qp.w - a.qm.w
        scl = convert(typeof(tmp), tmp / 2)
        Even(Quaternion(scl), Quaternion(-scl))
    else
        zero(a)
    end
end

function SimpleGA.project(a::Odd, n::Integer)
    return if (n == 1)
        convert(typeof(a), (a + a') / 2)
    elseif (n == 3)
        convert(typeof(a), (a - a') / 2)
    else
        zero(a)
    end
end

function LinearAlgebra.tr(a::Even)
    tmp = a.qp.w + a.qm.w
    return convert(typeof(tmp), tmp / 2)
end

function LinearAlgebra.dot(a::Even, b::Even)
    tmp = dot(a.qp, b.qp) + dot(a.qm, b.qm)
    return convert(typeof(tmp), tmp / 2)
end

function LinearAlgebra.dot(a::Odd, b::Odd)
    tmp = dot(a.qp, b.qm) + dot(a.qm, b.qp)
    return convert(typeof(tmp), tmp / 2)
end

#Exponentiation
SimpleGA.bivector_exp(a::Even) = Even(bivector_exp(a.qp), bivector_exp(a.qm))

Base.exp(a::Even) = Even(exp(a.qp), exp(a.qm))

#Comparison
function Base.isapprox(a::Even, b::Even; kwargs...)
    return isapprox(a.qp, b.qp; kwargs...) && isapprox(a.qm, b.qm; kwargs...)
end
function Base.isapprox(a::Odd, b::Odd; kwargs...)
    return isapprox(a.qp, b.qp; kwargs...) && isapprox(a.qm, b.qm; kwargs...)
end
Base.isequal(a::Even, b::Even) = isequal(a.qp, b.qp) && isequal(a.qm, b.qm)
Base.isequal(a::Odd, b::Odd) = isequal(a.qp, b.qp) && isequal(a.qm, b.qm)
