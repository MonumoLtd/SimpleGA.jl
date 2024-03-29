#=
Core code for the implementation of GA(3,0).
Underlying representation is with quaternions, though in this case we do not use the quaternion code as that would put in an unnecesary layer of indirection.
Instead we just use the hard-coded version of quaternions multiplication.
=#

struct Even{T<:Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

struct Odd{T<:Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

Even(w::Real, x::Real, y::Real, z::Real) = Even(promote(w, x, y, z)...)
Odd(w::Real, x::Real, y::Real, z::Real) = Odd(promote(w, x, y, z)...)

function Base.convert(::Type{Even{T}}, a::Even) where {T<:Real}
    return Even{T}(convert(T, a.w), convert(T, a.x), convert(T, a.y), convert(T, a.z))
end
function Base.convert(::Type{Odd{T}}, a::Odd) where {T<:Real}
    return Odd{T}(convert(T, a.w), convert(T, a.x), convert(T, a.y), convert(T, a.z))
end
function Base.promote_rule(::Type{Even{S}}, ::Type{Even{T}}) where {S<:Real,T<:Real}
    return Even{promote_type(S, T)}
end
function Base.promote_rule(::Type{Odd{S}}, ::Type{Odd{T}}) where {S<:Real,T<:Real}
    return Odd{promote_type(S, T)}
end
Base.zero(::Type{Even{T}}) where {T} = Even(zero(T), zero(T), zero(T), zero(T))
Base.zero(::Type{Odd{T}}) where {T} = Odd(zero(T), zero(T), zero(T), zero(T))
Base.one(::Type{Even{T}}) where {T} = Even(one(T), zero(T), zero(T), zero(T))
Base.zero(a::Even) = zero(typeof(a))
Base.zero(a::Odd) = zero(typeof(a))
Base.one(a::Even) = one(typeof(a))

# Addition / subtraction
Base.:(-)(a::Even) = Even(-a.w, -a.x, -a.y, -a.z)
Base.:(-)(a::Odd) = Odd(-a.w, -a.x, -a.y, -a.z)
Base.:(+)(a::Even, b::Even) = Even(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
Base.:(+)(a::Odd, b::Odd) = Odd(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
Base.:(-)(a::Even, b::Even) = Even(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
Base.:(-)(a::Odd, b::Odd) = Odd(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)

# Scalar addition / subtraction. Other cases are in GAcommon
# Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Real, a::Even) = Even(a.w + num, a.x, a.y, a.z)
Base.:(-)(num::Real, a::Even) = Even(-a.w + num, -a.x, -a.y, -a.z)

# Multiplication
Base.:(*)(num::Real, a::Even) = Even(num * a.w, num * a.x, num * a.y, num * a.z)
Base.:(*)(num::Real, a::Odd) = Odd(num * a.w, num * a.x, num * a.y, num * a.z)

function Base.:(*)(a::Even, b::Even)
    return Even(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

function Base.:(*)(a::Even, b::Odd)
    return Odd(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

function Base.:(*)(a::Odd, b::Even)
    return Odd(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

function Base.:(*)(a::Odd, b::Odd)
    return Even(
        -a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z,
        -a.w * b.x - a.x * b.w - a.y * b.z + a.z * b.y,
        -a.w * b.y - a.y * b.w - a.z * b.x + a.x * b.z,
        -a.w * b.z - a.z * b.w - a.x * b.y + a.y * b.x,
    )
end

# Reverse
LinearAlgebra.adjoint(a::Even) = Even(a.w, -a.x, -a.y, -a.z)
LinearAlgebra.adjoint(a::Odd) = Odd(-a.w, a.x, a.y, a.z)

# Grade and projection
function SimpleGA.project(a::Even, n::Integer)
    return if (n == 0)
        Even(a.w, zero(a.w), zero(a.w), zero(a.w))
    elseif (n == 2)
        Even(zero(a.w), a.x, a.y, a.z)
    else
        zero(a)
    end
end

function SimpleGA.project(a::Odd, n::Integer)
    return if (n == 3)
        Odd(a.w, zero(a.w), zero(a.w), zero(a.w))
    elseif (n == 1)
        Odd(zero(a.w), a.x, a.y, a.z)
    else
        zero(a)
    end
end

LinearAlgebra.tr(a::Even) = a.w
LinearAlgebra.dot(a::Even, b::Even) = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
LinearAlgebra.dot(a::Odd, b::Odd) = -a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z
LinearAlgebra.norm(a::Even) = sqrt(abs(dot(a, a)))
LinearAlgebra.norm(a::Odd) = sqrt(abs(dot(a, a)))

# Exponentiation
function SimpleGA.bivector_exp(a::Even)
    a = project(a, 2)
    nrm = sqrt(dot(a, -a))
    return if iszero(nrm)
        Even(one(a.w), zero(a.w), zero(a.w), zero(a.w))
    else
        cos(nrm) + sin(nrm) * a / nrm
    end
end

function Base.exp(a::Even)
    R = bivector_exp(a)
    return iszero(a.w) ? R : exp(a.w) * R
end

# Comparison.
function Base.isapprox(a::Even, b::Even; kwargs...)
    return isapprox(a.w, b.w; kwargs...) &&
           isapprox(a.x, b.x; kwargs...) &&
           isapprox(a.y, b.y; kwargs...) &&
           isapprox(a.z, b.z; kwargs...)
end
function Base.isapprox(a::Odd, b::Odd; kwargs...)
    return isapprox(a.w, b.w; kwargs...) &&
           isapprox(a.x, b.x; kwargs...) &&
           isapprox(a.y, b.y; kwargs...) &&
           isapprox(a.z, b.z; kwargs...)
end
function Base.isequal(a::Even, b::Even)
    return isequal(a.w, b.w) && isequal(a.x, b.x) && isequal(a.y, b.y) && isequal(a.z, b.z)
end
function Base.isequal(a::Odd, b::Odd)
    return isequal(a.w, b.w) && isequal(a.x, b.x) && isequal(a.y, b.y) && isequal(a.z, b.z)
end
