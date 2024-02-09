#=
Core code for the implementation of GA(3,1).
Underlying representation is with 2x2 complex matrices, though no matrices are formally constructed. The multiplication rules are hard-coded for efficiency.
Makes use of Julia's internal format.
=#

struct Even{T<:Real} <: Number
    c1::Complex{T}
    c2::Complex{T}
    c3::Complex{T}
    c4::Complex{T}
end

struct Odd{T<:Real} <: Number
    c1::Complex{T}
    c2::Complex{T}
    c3::Complex{T}
    c4::Complex{T}
end

Even(c1::Complex, c2::Complex, c3::Complex, c4::Complex) = Even(promote(c1, c2, c3, c4)...)
Odd(c1::Complex, c2::Complex, c3::Complex, c4::Complex) = Odd(promote(c1, c2, c3, c4)...)

function Base.convert(::Type{Even{T}}, a::Even) where {T<:Real}
    return Even{T}(
        convert(Complex{T}, a.c1),
        convert(Complex{T}, a.c2),
        convert(Complex{T}, a.c3),
        convert(Complex{T}, a.c4),
    )
end

function Base.convert(::Type{Odd{T}}, a::Odd) where {T<:Real}
    return Odd{T}(
        convert(Complex{T}, a.c1),
        convert(Complex{T}, a.c2),
        convert(Complex{T}, a.c3),
        convert(Complex{T}, a.c4),
    )
end
function Base.promote_rule(::Type{Even{S}}, ::Type{Even{T}}) where {S<:Real,T<:Real}
    return Even{promote_type(S, T)}
end
function Base.promote_rule(::Type{Odd{S}}, ::Type{Odd{T}}) where {S<:Real,T<:Real}
    return Odd{promote_type(S, T)}
end

function Base.zero(::Type{Even{T}}) where {T}
    z = zero(Complex{T})
    return Even(z, z, z, z)
end
function Base.zero(::Type{Odd{T}}) where {T}
    z = zero(Complex{T})
    return Odd(z, z, z, z)
end
function Base.one(::Type{Even{T}}) where {T}
    o = one(Complex{T})
    z = zero(Complex{T})
    return Even(o, z, z, o)
end
Base.zero(a::Even) = zero(typeof(a))
Base.zero(a::Odd) = zero(typeof(a))
Base.one(a::Even) = one(typeof(a))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.c1, -a.c2, -a.c3, -a.c4)
Base.:(-)(a::Odd) = Odd(-a.c1, -a.c2, -a.c3, -a.c4)
Base.:(+)(a::Even, b::Even) = Even(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(+)(a::Odd, b::Odd) = Odd(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(-)(a::Even, b::Even) = Even(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)
Base.:(-)(a::Odd, b::Odd) = Odd(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)

#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Real, a::Even) = Even(a.c1 + num, a.c2, a.c3, a.c4 + num)
Base.:(-)(num::Real, a::Even) = Even(-a.c1 + num, -a.c2, -a.c3, -a.c4 + num)

#Multiplication
Base.:(*)(num::Real, a::Even) =Even(num * a.c1, num * a.c2, num * a.c3, num * a.c4)
Base.:(*)(num::Real, a::Odd) = Odd(num * a.c1, num * a.c2, num * a.c3, num * a.c4)

function Base.:(*)(a::Even, b::Even)
    return Even(
        a.c1 * b.c1 + a.c2 * b.c3,
        a.c1 * b.c2 + a.c2 * b.c4,
        a.c3 * b.c1 + a.c4 * b.c3,
        a.c4 * b.c4 + a.c3 * b.c2,
    )
end

function Base.:(*)(a::Even, b::Odd)
    return Odd(
        a.c1 * b.c1 + a.c2 * b.c3,
        a.c1 * b.c2 + a.c2 * b.c4,
        a.c3 * b.c1 + a.c4 * b.c3,
        a.c4 * b.c4 + a.c3 * b.c2,
    )
end

function Base.:(*)(a::Odd, b::Even)
    return Odd(
        a.c1 * conj(b.c4) - a.c2 * conj(b.c2),
        -a.c1 * conj(b.c3) + a.c2 * conj(b.c1),
        a.c3 * conj(b.c4) - a.c4 * conj(b.c2),
        a.c4 * conj(b.c1) - a.c3 * conj(b.c3),
    )
end

function Base.:(*)(a::Odd, b::Odd)
    return Even(
        -a.c1 * conj(b.c4) + a.c2 * conj(b.c2),
        a.c1 * conj(b.c3) - a.c2 * conj(b.c1),
        -a.c3 * conj(b.c4) + a.c4 * conj(b.c2),
        -a.c4 * conj(b.c1) + a.c3 * conj(b.c3),
    )
end

#Reverse
LinearAlgebra.adjoint(a::Even) = Even(a.c4, -a.c2, -a.c3, a.c1)
LinearAlgebra.adjoint(a::Odd) = Odd(conj(a.c1), conj(a.c3), conj(a.c2), conj(a.c4))

#Grade and projection
function SimpleGA.project(a::Even, n::Integer)
    tmp = (a.c1 + a.c4)
    tra = convert(typeof(tmp), tmp / 2)
    return if (n == 0)
        real(tra) * one(a)
    elseif (n == 2)
        convert(typeof(a), (a - a') / 2)
    elseif (n == 4)
        Even(imag(tra) * im, zero(a.c1), zero(a.c1), imag(tra) * im)
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
        return zero(a)
    end
end

function LinearAlgebra.tr(a::Even)
    tmp = (real(a.c1 + a.c4))
    return convert(typeof(tmp), tmp / 2)
end

function LinearAlgebra.dot(a::Even, b::Even)
    tmp = real(a.c1 * b.c1 + a.c2 * b.c3 + a.c4 * b.c4 + a.c3 * b.c2)
    return convert(typeof(tmp), tmp / 2)
end
function LinearAlgebra.dot(a::Odd, b::Odd)
    tmp = real(
        -a.c1 * conj(b.c4) + a.c2 * conj(b.c2) - a.c4 * conj(b.c1) + a.c3 * conj(b.c3)
    )
    return convert(typeof(tmp), tmp / 2)
end
LinearAlgebra.norm(a::Even) = sqrt(abs(dot(a, a)))
LinearAlgebra.norm(a::Odd) = sqrt(abs(dot(a, a)))

"""
    _as_even(x::Complex) -> Even

Wrap a complex number `x` as an Even type.

We don't expose this so as to avoid leaking our representation abstraction in the definition
of +, - and *.
"""
_as_even(x::Complex) = Even(x, zero(x), zero(x), x)

#Exponentiation
function SimpleGA.bivector_exp(a::Even)
    a = project(a, 2)
    aa = a * a
    fct = sqrt((aa.c1 + aa.c4) / 2)
    return if iszero(fct)
        1 + a
    else
        _as_even(cosh(fct)) + _as_even(sinh(fct) / fct) * a
    end
end

function Base.exp(a::Even)
    R = bivector_exp(a)
    fct = (a.c1 + a.c4) / 2
    return if iszero(fct)
        R
    else
        _as_even(exp(fct)) * R
    end
end

# Comparison.
function Base.isapprox(a::Even, b::Even; kwargs...)
    return isapprox(a.c1, b.c1; kwargs...) &&
           isapprox(a.c2, b.c2; kwargs...) &&
           isapprox(a.c3, b.c3; kwargs...) &&
           isapprox(a.c4, b.c4; kwargs...)
end
function Base.isapprox(a::Odd, b::Odd; kwargs...)
    return isapprox(a.c1, b.c1; kwargs...) &&
           isapprox(a.c2, b.c2; kwargs...) &&
           isapprox(a.c3, b.c3; kwargs...) &&
           isapprox(a.c4, b.c4; kwargs...)
end
function Base.isequal(a::Even, b::Even)
    return isequal(a.c1, b.c1) &&
           isequal(a.c2, b.c2) &&
           isequal(a.c3, b.c3) &&
           isequal(a.c4, b.c4)
end
function Base.isequal(a::Odd, b::Odd)
    return isequal(a.c1, b.c1) &&
           isequal(a.c2, b.c2) &&
           isequal(a.c3, b.c3) &&
           isequal(a.c4, b.c4)
end
