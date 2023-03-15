module Quaternions
#=
Quaternion code. This is called by some of the later GA implementations. 
The core mirrors much of the GA code structure.
For completeness we have defined a division operation for quaternions as they are a division algebra.
=#

using LinearAlgebra

import ..project
import ..expb


export real_part, imag_part, expb, Quaternion

struct Quaternion{T<:Real}
    w::T
    x::T
    y::T
    z::T
end

function Base.zero(::Type{Quaternion{T}}) where {T}
    return Quaternion{T}(zero(T), zero(T), zero(T), zero(T))
end

Base.:(-)(a::Quaternion) = Quaternion(-a.w, -a.x, -a.y, -a.z)
Base.:(+)(a::Quaternion, b::Quaternion) = Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
Base.:(-)(a::Quaternion, b::Quaternion) = Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
Base.:(+)(x::Number, a::Quaternion) = Quaternion(x + a.w, a.x, a.y, a.z)
Base.:(-)(x::Number, a::Quaternion) = Quaternion(x - a.w, -a.x, -a.y, -a.z)
Base.:(+)(a::Quaternion, x::Real) = Quaternion(a.w + x, a.x, a.y, a.z)
Base.:(-)(a::Quaternion, x::Real) = Quaternion(a.w - x, a.x, a.y, a.z)

#Multiplication
Base.:(*)(x::Real, a::Quaternion) = Quaternion(x * a.w, x * a.x, x * a.y, x * a.z)
Base.:(*)(a::Quaternion, x::Real) = x*a

function Base.:(*)(a::Quaternion, b::Quaternion)
    return Quaternion(
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
        a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
        a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
        a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
    )
end

Base.:(/)(a::Quaternion, x::Real) = (1/x)*a
Base.:(/)(a::Quaternion, b::Quaternion) = a*conj(b) / dot(b,conj(b))
LinearAlgebra.norm(a::Quaternion) = sqrt(a.w^2+a.x^2+a.y^2+a.z^2)

#Reverse
Base.conj(a::Quaternion) = Quaternion(a.w, -a.x, -a.y, -a.z)


#Projection operations
LinearAlgebra.dot(a::Quaternion, b::Quaternion) = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
Base.real(a::Quaternion) = a.w
real_part(a::Quaternion) = Quaternion(a.w,zero(a.w),zero(a.w),zero(a.w))
imag_part(a::Quaternion) = Quaternion(zero(a.w),a.x,a.y,a.z)

#Exponentiation
function expb(a::Quaternion)
    a = imag_part(a)
    nrm = norm(a)
    if iszero(nrm)
        return Quaternion(one(nrm),0,0,0)
    else
        return cos(nrm) + sin(nrm)*a/nrm
    end
end

function exp(a::Quaternion)
    R = expb(a)
    if iszero(a.w)
        return R
    else 
        return exp(a.w)*R
    end
end



#Additional Functions
Base.isapprox(a::Quaternion, b::Quaternion) = isapprox(a.w,b.w) && isapprox(a.x,b.x) && isapprox(a.y,b.y) && isapprox(a.z,b.z) 

function Base.isapprox(a::Quaternion{T}, b::Quaternion{T}; kwargs...) where {T}
    return (
        isapprox(a.w, b.w; kwargs...) &&
        isapprox(a.x, b.x; kwargs...) &&
        isapprox(a.y, b.y; kwargs...) &&
        isapprox(a.z, b.z; kwargs...)
    )
end


function Base.show(io::IO, ::MIME"text/plain", a::Quaternion)
    print(io, "", string(a.w) * " + " * string(a.x) * "i + " * string(a.y) * "j + " * string(a.z) * "k")
end


function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Quaternion})
    n= length(mvs)
    println(io,n,"-element Vector{Quaternion}")
    for i in eachindex(mvs)
    println(io, " ", mvs[i])
    end
end

end #Module