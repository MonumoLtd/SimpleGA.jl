"""
Core code for the implementation of GA(2,0)

Underlying representation is with complex numbers, so essentially a wrapper for Julia's
internal ComplexF{T} formats.
"""

struct Even{T<:Real} <: Number
    c1::Complex{T}
end

struct Odd{T<:Real} <: Number
    c1::Complex{T}
end

Base.convert(::Type{Even{T}}, a::Even) where {T<:Real} = Even{T}(convert(Complex{T}, a.c1))
Base.convert(::Type{Odd{T}}, a::Odd) where {T<:Real} = Odd{T}(convert(Complex{T}, a.c1))
function Base.promote_rule(::Type{Even{S}}, ::Type{Even{T}}) where {S<:Real,T<:Real}
    return Even{promote_rule(S, T)}
end
function Base.promote_rule(::Type{Odd{S}}, ::Type{Odd{T}}) where {S<:Real,T<:Real}
    return Odd{promote_rule(S, T)}
end
Base.zero(a::Even) = Even(zero(a.c1))
Base.zero(a::Odd) = Odd(zero(a.c1))
Base.one(a::Even) = Even(one(a.c1))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.c1)
Base.:(-)(a::Odd) = Odd(-a.c1)
Base.:(+)(a::Even, b::Even) = Even(a.c1 + b.c1)
Base.:(+)(a::Odd, b::Odd) = Odd(a.c1 + b.c1)
Base.:(-)(a::Even, b::Even) = Even(a.c1 - b.c1)
Base.:(-)(a::Odd, b::Odd) = Odd(a.c1 - b.c1)

#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number, a::Even) = Even(a.c1 + num)
Base.:(-)(num::Number, a::Even) = Even(-a.c1 + num)

#Multiplication
Base.:(*)(num::Number, a::Even) = Even(num * a.c1)
Base.:(*)(num::Number, a::Odd) = Odd(num * a.c1)
Base.:(*)(a::Even, b::Even) = Even(a.c1 * b.c1)
Base.:(*)(a::Even, b::Odd) = Odd(conj(a.c1) * b.c1)
Base.:(*)(a::Odd, b::Even) = Odd(a.c1 * b.c1)
Base.:(*)(a::Odd, b::Odd) = Even(conj(a.c1) * b.c1)

#Allow this here as complex numbers are a division algebra.
Base.:(/)(a::Even, b::Even) = Even(a.c1 / b.c1)

#Reverse
LinearAlgebra.adjoint(a::Even) = Even(conj(a.c1))
LinearAlgebra.adjoint(a::Odd) = a

#Grade and projection
function SimpleGA.project(a::Even, n::Integer)
    return if (n == 0)
        real(a.c1) * one(a)
    elseif (n == 2)
        Even(imag(a.c1) * im)
    else
        zero(a)
    end
end

SimpleGA.project(a::Odd, n::Integer) = n == 1 ? a : zero(a)

LinearAlgebra.tr(a::Even) = real(a.c1)
LinearAlgebra.dot(a::Even, b::Even) = real(a.c1 * b.c1)
LinearAlgebra.dot(a::Odd, b::Odd) = real(conj(a.c1) * b.c1)

#Exponentiation
Base.exp(a::Even) = Even(exp(a.c1))

SimpleGA.bivector_exp(a::Even) = Even(exp(im * a.c1.im))

# Comparison
Base.isapprox(a::Even, b::Even; kwargs...) = isapprox(a.c1, b.c1; kwargs...)
Base.isapprox(a::Odd, b::Odd; kwargs...) = isapprox(a.c1, b.c1; kwargs...)
Base.isequal(a::Even, b::Even) = isequal(a.c1, b.c1)
Base.isequal(a::Odd, b::Odd) = isequal(a.c1, b.c1)
