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

Base.zero(a::Even) = Even(zero(a.c1), zero(a.c2), zero(a.c3), zero(a.c4))
Base.zero(a::Odd) = Odd(zero(a.c1), zero(a.c2), zero(a.c3), zero(a.c4))
Base.one(a::Even) = Even(one(a.c1), zero(a.c2), zero(a.c3), one(a.c4))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.c1, -a.c2, -a.c3, -a.c4)
Base.:(-)(a::Odd) = Odd(-a.c1, -a.c2, -a.c3, -a.c4)
Base.:(+)(a::Even, b::Even) = Even(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(+)(a::Odd, b::Odd) = Odd(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(-)(a::Even, b::Even) = Even(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)
Base.:(-)(a::Odd, b::Odd) = Odd(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)

#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number, a::Even) = Even(a.c1 + num, a.c2, a.c3, a.c4 + num)
Base.:(-)(num::Number, a::Even) = (-a.c1 + num, -a.c2, -a.c3, -a.c4 + num)

#Multiplication
Base.:(*)(num::Number, a::Even) = Even(num * a.c1, num * a.c2, num * a.c3, num * a.c4)
Base.:(*)(num::Number, a::Odd) = Odd(num * a.c1, num * a.c2, num * a.c3, num * a.c4)

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
       - a.c3 * conj(b.c4) + a.c4 * conj(b.c2),
        -a.c4 * conj(b.c1) + a.c3 * conj(b.c3),
    )
end

#Reverse
LinearAlgebra.adjoint(a::Even) = Even(a.c4, -a.c2, -a.c3, a.c1)
LinearAlgebra.adjoint(a::Odd) = Odd(conj(a.c1), conj(a.c3), conj(a.c2), conj(a.c4))

#Grade and projection
function GeometricAlgebra.project(a::Even, n::Integer)
    tra = (a.c1 + a.c4) / 2
    return if (n == 0)
        real(tra) * one(a)
    elseif (n == 2)
        convert(typeof(a),(a - a') / 2)
    elseif (n == 4)
        Even((tra - conj(tra)) / 2, zero(a.c1), zero(a.c1), (tra - conj(tra)) / 2)
    else
        zero(a)
    end
end

function GeometricAlgebra.project(a::Odd, n::Integer)
    return if (n == 1)
        convert(typeof(a),(a + a') / 2)
    elseif (n == 3)
        convert(typeof(a),(a - a') / 2)
    else
        return zero(a)
    end
end

LinearAlgebra.tr(a::Even) = (real(a.c1 + a.c4)) / 2

function LinearAlgebra.dot(a::Even, b::Even)
    return real(a.c1 * b.c1 + a.c2 * b.c3 + a.c4 * b.c4 + a.c3 * b.c2) / 2
end
function LinearAlgebra.dot(a::Odd, b::Odd)
    return real(
        -a.c1 * conj(b.c4) + a.c2 * conj(b.c2) - a.c4 * conj(b.c1) + a.c3 * conj(b.c3)
    ) / 2
end

#Exponentiation
function GeometricAlgebra.bivector_exp(a::Even)
    a = project(a, 2)
    aa = a * a
    fct = sqrt((aa.c1 + aa.c4) / 2)
    return if iszero(fct)
        1 + a
    else
        cosh(fct) + sinh(fct) / fct * a
    end
end

function Base.exp(a::Even)
    R = bivector_exp(a)
    fct = (a.c1 + a.c4) / 2
    return if iszero(fct)
        R
    else
        exp(fct) * R
    end
end

#Comparison. Uses default tolerances
function Base.isapprox(a::Even, b::Even)
    return isapprox(a.c1, b.c1) &&
           isapprox(a.c2, b.c2) &&
           isapprox(a.c3, b.c3) &&
           isapprox(a.c4, b.c4)
end
function Base.isapprox(a::Odd, b::Odd)
    return isapprox(a.c1, b.c1) &&
           isapprox(a.c2, b.c2) &&
           isapprox(a.c3, b.c3) &&
           isapprox(a.c4, b.c4)
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
