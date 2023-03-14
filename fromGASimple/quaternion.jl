#=
Quaternion code. This is called by some of the later GA implementations. 
The core mirrors much of the GA code structure.
For completeness we have defined a division operation for quaternions as they are a division algebra.
Pretty-typing has not been implented here.
=#


import Base.:*
import Base.:+
import Base.:-
import Base.:/
import Base.exp
import Base.real
import Base.conj
import LinearAlgebra.norm
import LinearAlgebra.dot
import Base.show


struct Quaternion
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

function quat(w,x,y,z)
    Quaternion(w,x,y,z)
end

const qzero = Quaternion(0,0,0,0)

#Maths operations
function -(a::Quaternion)
    Quaternion(-a.w,-a.x,-a.y,-a.z)
end

function +(a::Quaternion, b::Quaternion)
    Quaternion(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function -(a::Quaternion,b::Quaternion)
    Quaternion(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

function +(num::Number,a::Quaternion)
    Quaternion(a.w+convert(Float64,num), a.x, a.y, a.z)
end

function -(num::Number,a::Quaternion)
    Quaternion(-a.w+convert(Float64,num), -a.x, -a.y, -a.z)
end

function +(a::Quaternion,num::Number)
    Quaternion(a.w+convert(Float64,num), a.x, a.y, a.z)
end

function -(a::Quaternion,num::Number)
    Quaternion(a.w-convert(Float64,num), a.x, a.y, a.z)
end

function *(num::Number,a::Quaternion)
    num = convert(Float64,num)
    Quaternion(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(a::Quaternion,num::Number)
    num = convert(Float64,num)
    Quaternion(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(a::Quaternion, b::Quaternion)
    Quaternion(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
             a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
             a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
             a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function /(a::Quaternion,num::Number)
    num = convert(Float64,1/num)
    Quaternion(num*a.w, num*a.x, num*a.y, num*a.z)
end

function /(a::Quaternion,b::Quaternion)
    a*conj(b) / dot(b,conj(b))
end

function norm(a::Quaternion)
    return sqrt(a.w^2+a.x^2+a.y^2+a.z^2)
end

#Reverse
function conj(a::Quaternion)
    Quaternion(a.w,-a.x, -a.y, -a.z)
end

#Projection operations
function dot(a::Quaternion, b::Quaternion)
    a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
end

function real(a::Quaternion)
    a.w
end

function real_part(a::Quaternion)
    Quaternion(a.w,0,0,0)
end

function imag_part(a::Quaternion)
    Quaternion(0,a.x,a.y,a.z)
end


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

function Base.show(io::IO, a::Quaternion)
    print(string(a.w) * " + " * string(a.x) * "i + " * string(a.y) * "j + " * string(a.z) * "k")
end


function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Quaternion})
    n= length(mvs)
    println(io,n,"-element Vector{Quaternion}")
    for i in eachindex(mvs)
    println(io, " ", mvs[i])
    end
end