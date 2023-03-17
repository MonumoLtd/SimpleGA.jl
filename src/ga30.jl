"""
    GA30

Module representing the GA(3, 0) geometric algebra.
Underlying rep is quaternions, but not used explicitly here to keep this code self-contained.
"""

module GA30

using LinearAlgebra

include("core30.jl")
include("common.jl")

#Basis
const e1 = Odd{Int8}(0,1,0,0)
const e2 = Odd{Int8}(0,0,1,0)
const e3 = Odd{Int8}(0,0,0,1)
const I3 = Odd{Int8}(1,0,0,0)

bas30 = [e1,e2,e3]
export bas30


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)


function mvtype(a::Even)
    res=""
    tp = approxzero(a.w) ? "" : " + " * string(a.w)
    res *= tp
    tp = approxzero(a.x) ? "" : " + " * string(-a.x) * "e2e3"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(-a.y) * "e3e1"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(-a.z) * "e1e2"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

function mvtype(a::Odd)
    res=""
    tp = approxzero(a.x) ? "" : " + " * string(a.x) * "e1"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(a.y) * "e2"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(a.z) * "e3"
    res *= tp
    tp = approxzero(a.w) ? "" : " + " * string(a.w) * "e123"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

include("show.jl")

end

#=
struct MV{P,T} <: MultiVector{P,T}
    q::Quaternion{T}
end
MV{P}(q::Quaternion{T}) where {P,T} = MV{P,T}(q)
function MV{P,T}(w, x, y, z) where {P<:Parity,T<:Real}
    return MV{P,typeof(w)}(Quaternion{T}(w, x, y, z))
end
function MV{P}(w, x, y, z) where {P<:Parity}
    w, x, y, z = promote(w, x, y, z)
    return MV{P,typeof(w)}(w, x, y, z)
end

# Construction
Base.zero(::MV{P,T}) where {P<:Parity,T} = MV{P}(zero(Quaternion{T}))

# Comparison
Base.isapprox(a::MV{P}, b::MV{P}; kwargs...) where {P} = isapprox(a.q, b.q; kwargs...)

function Base.show(io::IO, a::MV)
    parts = if iseven(a)
        vcat(
            string_parts(a.q.w),
            string_parts(-a.q.x, "e2e3"),
            string_parts(-a.q.y, "e3e1"),
            string_parts(-a.q.z, "e1e2"),
        )
    else  # isodd(a)
        vcat(
            string_parts(a.q.x, "e1"),
            string_parts(a.q.y, "e2"),
            string_parts(a.q.z, "e3"),
            string_parts(a.q.w, "I3"),
        )
    end
    return print(io, join(parts, " + "))
end

# Unary negation.
Base.:(-)(a::MV{P}) where {P} = MV{P}(-a.q)

# Addition & subtraction only defined between matching parities.
Base.:(+)(a::MV{P}, b::MV{P}) where {P} = MV{P}(a.q + b.q)
Base.:(-)(a::MV{P}, b::MV{P}) where {P} = MV{P}(a.q - b.q)

# Multiplication rules depend on parity.
Base.:(*)(a::MV{Even}, b::MV{Even}) = MV{Even}(a.q * b.q)
Base.:(*)(a::MV{Even}, b::MV{Odd}) = MV{Odd}(a.q * b.q)
Base.:(*)(a::MV{Odd}, b::MV{Even}) = MV{Odd}(a.q * b.q)
Base.:(*)(a::MV{Odd}, b::MV{Odd}) = MV{Even}(-a.q * b.q)

# Addition and subtraction with scalars.
Base.:(+)(x::Real, a::MV{Even}) = MV{Even}(x + a.q)
Base.:(+)(a::MV{Even}, x::Real) = MV{Even}(a.q + x)
Base.:(-)(x::Real, a::MV{Even}) = MV{Even}(x - a.q)
Base.:(-)(a::MV{Even}, x::Real) = MV{Even}(a.q - x)

# Multiplication with scalars.
Base.:(*)(x::Real, a::MV{P}) where {P} = MV{P}(x * a.q)
Base.:(*)(a::MV{P}, x::Real) where {P} = MV{P}(a.q * x)

# TODO division

Base.reverse(a::MV{Even}) = MV{Even}(reverse(a.q))
Base.reverse(a::MV{Odd}) = MV{Odd}(-reverse(a.q))

function GADraft.project(a::MV{Even,T}, grade::Integer) where {T}
    return if grade == 0
        MV{Even,T}(a.q.w, zero(T), zero(T), zero(T))
    elseif grade == 2
        MV{Even,T}(zero(T), a.q.x, a.q.y, a.q.z)
    else
        # Any other grade won't have a contribution.
        # TODO should we throw exceptions for non-sensical grades?
        #   (e.g. negative or larger than 2)
        zero(a)
    end
end

function GADraft.project(a::MV{Odd,T}, grade::Integer) where {T}
    return if grade == 3
        MV{Odd,T}(a.q.w, zero(T), zero(T), zero(T))
    elseif grade == 1
        MV{Odd,T}(zero(T), a.q.x, a.q.y, a.q.z)
    else
        # Any other grade won't have a contribution.
        # TODO should we throw exceptions for non-sensical grades?
        zero(a)
    end
end

GADraft.scalar(a::MV{Even}) = scalar(a.q)
GADraft.scalar(a::MV{Even}, b::MV{Even}) = scalar(a.q, b.q)
GADraft.scalar(a::MV{Odd}, b::MV{Odd}) = -scalar(a.q, b.q)

# Default basis uses Float64.
const e1 = MV{Odd,Float64}(0.0, 1.0, 0.0, 0.0)
const e2 = MV{Odd,Float64}(0.0, 0.0, 1.0, 0.0)
const e3 = MV{Odd,Float64}(0.0, 0.0, 0.0, 1.0)
const I3 = MV{Odd,Float64}(1.0, 0.0, 0.0, 0.0)

end  # module GA30
=#