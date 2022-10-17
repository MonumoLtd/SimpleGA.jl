struct Quaternion{T<:Real}
    w::T
    x::T
    y::T
    z::T
end

function Base.zero(::Type{Quaternion{T}}) where {T}
    return Quaternion{T}(zero(T), zero(T), zero(T), zero(T))
end

# Comparison
function Base.isapprox(a::Quaternion{T}, b::Quaternion{T}; kwargs...) where {T}
    return (
        isapprox(a.w, b.w; kwargs...) &&
        isapprox(a.x, b.x; kwargs...) &&
        isapprox(a.y, b.y; kwargs...) &&
        isapprox(a.z, b.z; kwargs...)
    )
end

Base.:(-)(a::Quaternion) = Quaternion(-a.w, -a.x, -a.y, -a.z)

function Base.:(+)(a::Quaternion{T}, b::Quaternion{T}) where {T}
    return Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function Base.:(-)(a::Quaternion{T}, b::Quaternion{T}) where {T}
    return Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

Base.:(+)(x::Real, a::Quaternion) = Quaternion(x + a.w, a.x, a.y, a.z)
Base.:(-)(x::Real, a::Quaternion) = Quaternion(x - a.w, -a.x, -a.y, -a.z)
Base.:(+)(a::Quaternion, x::Real) = Quaternion(a.w + x, a.x, a.y, a.z)
Base.:(-)(a::Quaternion, x::Real) = Quaternion(a.w - x, a.x, a.y, a.z)

function Base.:(*)(x::Real, a::Quaternion{T}) where {T}
    x = convert(T, x)
    return Quaternion(x * a.w, x * a.x, x * a.y, x * a.z)
end

function Base.:(*)(a::Quaternion{T}, x::Real) where {T}
    x = convert(T, x)
    return Quaternion(x * a.w, x * a.x, x * a.y, x * a.z)
end

function Base.:(*)(a::Quaternion, b::Quaternion)
    return Quaternion(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

function Base.:(/)(a::Quaternion{T}, x::Real) where {T}
    x = convert(T, inverse(x))
    return Quaternion(x * a.w, x * a.x, x * a.y, x * a.z)
end

function Base.:(/)(a::Quaternion, b::Quaternion)
    return a * reverse(b) / scalar(b, reverse(b))
end

Base.reverse(a::Quaternion) = Quaternion(a.w, -a.x, -a.y, -a.z)

#Projection operations
GADraft.scalar(a::Quaternion) = a.w

GADraft.scalar(a::Quaternion, b::Quaternion) = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z

function GADraft.project(a::Quaternion{T}, grade::Integer) where {T}
    return if (grade == 0)
        Quaternion{T}(a.w, zero(T), zero(T), zero(T))
    elseif (grade == 2)
        Quaternion{T}(zero(T), a.x, a.y, a.z)
    else
        zero(a)
    end
end
