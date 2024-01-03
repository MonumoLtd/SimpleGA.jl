#=
Core code for the implementation of GA(3,3).
Work using self-dual and anti-self-dual decomposition. Base element is a 4x4 matrix built on Static Arrays library.
Useful algebra for line geometry.
=#

struct Even{T<:Real} <: Number
    p::SMatrix{4,4,T,16}
    m::SMatrix{4,4,T,16}
end

struct Odd{T<:Real} <: Number
    p::SMatrix{4,4,T,16}
    m::SMatrix{4,4,T,16}
end

id4 = SMatrix{4,4,Int8}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
J4 = SMatrix{4,4,Int8}([0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0])

function Base.convert(::Type{Even{T}}, a::Even) where {T<:Real}
    return Even{T}(convert(SMatrix{4,4,T,16}, a.p), convert(SMatrix{4,4,T,16}, a.m))
end

function Base.convert(::Type{Odd{T}}, a::Odd) where {T<:Real}
    return Odd{T}(convert(SMatrix{4,4,T,16}, a.p), convert(SMatrix{4,4,T,16}, a.m))
end
function Base.promote_rule(::Type{Even{S}}, ::Type{Even{T}}) where {S<:Real,T<:Real}
    return Even{promote_rule(S, T)}
end
function Base.promote_rule(::Type{Odd{S}}, ::Type{Odd{T}}) where {S<:Real,T<:Real}
    return Odd{promote_rule(S, T)}
end

Base.zero(a::Even) = Even(zero(a.p), zero(a.m))
Base.zero(a::Odd) = Odd(zero(a.p), zero(a.m))
Base.one(a::Even) = Even(one(a.p), one(a.m))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.p, -a.m)
Base.:(-)(a::Odd) = Odd(-a.p, -a.m)
Base.:(+)(a::Even, b::Even) = Even(a.p + b.p, a.m + b.m)
Base.:(+)(a::Odd, b::Odd) = Odd(a.p + b.p, a.m + b.m)
Base.:(-)(a::Even, b::Even) = Even(a.p - b.p, a.m - b.m)
Base.:(-)(a::Odd, b::Odd) = Odd(a.p - b.p, a.m - b.m)

#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number, a::Even) = Even(a.p + num * one(a.p), a.m + num * one(a.m))
Base.:(-)(num::Number, a::Even) = Even(-a.p + num * one(a.p), -a.m + num * one(a.m))

#Multiplication
Base.:(*)(num::Number, a::Even) = Even(num * a.p, num * a.m)
Base.:(*)(num::Number, a::Odd) = Odd(num * a.p, num * a.m)
Base.:(*)(a::Even, b::Even) = Even(a.p * b.p, a.m * b.m)
Base.:(*)(a::Even, b::Odd) = Odd(a.p * b.p, a.m * b.m)
Base.:(*)(a::Odd, b::Even) = Odd(a.p * b.m, a.m * b.p)
Base.:(*)(a::Odd, b::Odd) = Even(a.p * b.m, a.m * b.p)

#Reverse
function LinearAlgebra.adjoint(a::Even)
    return Even(-J4 * transpose(a.m) * J4, -J4 * transpose(a.p) * J4)
end

function LinearAlgebra.adjoint(a::Odd)
    return Odd(-J4 * transpose(a.p) * J4, -J4 * transpose(a.m) * J4)
end

#Grade and projection
function SimpleGA.project(a::Even, n::Integer)
    return if (n == 0)
        scl = (tr(a.p) + tr(a.m)) / 8
        convert(typeof(a), scl * one(a))
    elseif (n == 2)
        scl = (tr(a.p) - tr(a.m)) / 8
        convert(typeof(a), (a - a') / 2 - scl * Even(one(a.p), -one(a.m)))
    elseif (n == 4)
        scl = (tr(a.p) + tr(a.m)) / 8
        convert(typeof(a), (a + a') / 2 - scl * one(a))
    elseif (n == 6)
        scl = (tr(a.p) - tr(a.m)) / 8
        convert(typeof(a), scl * Even(one(a.p), -one(a.m)))
    else
        zero(a)
    end
end

function SimpleGA.project(a::Odd, n::Integer)
    return if (n == 3)
        convert(typeof(a), (a - a') / 2)
    elseif (n == 1)
        tmp = (a + a') * e3 / 2
        convert(typeof(a), (project(tmp, 0) + project(tmp, 2)) * e3)
    elseif (n == 5)
        tmp = (a + a') * e3 / 2
        convert(typeof(a), (project(tmp, 4) + project(tmp, 6)) * e3)
    else
        zero(a)
    end
end

function LinearAlgebra.tr(a::Even)
    tmp = tr(a.p) + tr(a.m)
    return convert(typeof(tmp), tmp / 8)
end

function LinearAlgebra.dot(a::Even, b::Even)
    tmp = tr(a.p * b.p) + tr(a.m * b.m)
    return convert(typeof(tmp), tmp / 8)
end

function LinearAlgebra.dot(a::Odd, b::Odd)
    tmp = tr(a.p * b.m) + tr(a.m * b.p)
    return convert(typeof(tmp), tmp / 8)
end

#Exponentiation
Base.exp(a::Even) = Even(exp(a.p), exp(a.m))

#TODO Any improvement here?
function SimpleGA.bivector_exp(a::Even)
    a = project(a, 2)
    R = exp(a)
    delt = R * R' - 1
    return (1 - 0.5 * delt + 0.375 * delt * delt) * R
end

#Comparison
#StaticArrays does seem to lose some accuracy.
function Base.isapprox(a::Even, b::Even; kwargs...)
    return isapprox(a.p, b.p; kwargs...) && isapprox(a.m, b.m; kwargs...)
end
function Base.isapprox(a::Odd, b::Odd; kwargs...)
    return isapprox(a.p, b.p; kwargs...) && isapprox(a.m, b.m; kwargs...)
end
Base.isequal(a::Even, b::Even) = isequal(a.p, b.p) && isequal(a.m, b.m)
Base.isequal(a::Odd, b::Odd) = isequal(a.p, b.p) && isequal(a.m, b.m)
