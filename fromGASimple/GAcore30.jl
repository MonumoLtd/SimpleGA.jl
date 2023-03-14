#=
Core code for the implementation of GA(3,0).
Underlying representation is with quaternions, though in this case we do not use the quaternion code as that would put in an unnecesary layer of indirection.
Instead we just use the hard-coded version of quaternions multiplication.
=#

import Base.:*
import Base.:+
import Base.:-
import Base.:/
import Base.exp
import LinearAlgebra.tr
import LinearAlgebra.dot
import LinearAlgebra.adjoint
import ..project
import ..expb

struct MVeven
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

struct MVodd
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

#Addition / subtraction
function -(a::MVeven)
    MVeven(-a.w,-a.x,-a.y,-a.z)
end

function -(a::MVodd)
    MVodd(-a.w,-a.x,-a.y,-a.z)
end

function +(a::MVeven, b::MVeven)
    MVeven(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function (-)(a::MVeven,b::MVeven)
    MVeven(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

function +(a::MVodd, b::MVodd)
    MVodd(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function -(a::MVodd, b::MVodd)
    MVodd(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::MVeven)
    MVeven(a.w+convert(Float64,num), a.x, a.y, a.z)
end

function -(num::Number,a::MVeven)
    MVeven(-a.w+convert(Float64,num), -a.x, -a.y, -a.z)
end

#Multiplication
function *(num::Number,a::MVeven)
    num = convert(Float64,num)
    MVeven(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(num::Number,a::MVodd)
    num = convert(Float64,num)
    MVodd(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(a::MVeven, b::MVeven)
    MVeven(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
             a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
             a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
             a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::MVeven, b::MVodd)
    MVodd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::MVodd, b::MVeven)
    MVodd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::MVodd, b::MVodd)
    MVeven(-a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z,
            -a.w*b.x - a.x*b.w - a.y*b.z + a.z*b.y,
            -a.w*b.y - a.y*b.w - a.z*b.x + a.x*b.z,
            -a.w*b.z - a.z*b.w - a.x*b.y + a.y*b.x )
end

#Division by a real
function /(a::MVeven,num::Number)
    num = convert(Float64,1/num)
    MVeven(num*a.w, num*a.x, num*a.y, num*a.z)
end

function /(a::MVodd,num::Number)
    num = convert(Float64,1/num)
    MVodd(num*a.w, num*a.x, num*a.y, num*a.z)
end

#Reverse
function adjoint(a::MVeven)
    MVeven(a.w,-a.x, -a.y, -a.z)
end

function adjoint(a::MVodd)
    MVodd(-a.w,a.x, a.y, a.z)
end

#Grade and projection
function project(a::MVeven,n::Integer)
    if (n == 0)
        return MVeven(a.w,0,0,0)
    elseif (n == 2)
        return MVeven(0,a.x,a.y,a.z)
    else
        return MVeven(0,0,0,0)
    end
end

function project(a::MVodd,n::Integer)
    if (n == 3)
        return MVodd(a.w,0,0,0)
    elseif (n == 1)
        return MVodd(0,a.x,a.y,a.z)
    else
        return MVodd(0,0,0,0)
    end
end

function tr(a::MVeven)
    a.w
end

function dot(a::MVeven, b::MVeven)
    a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
end

function dot(a::MVodd, b::MVodd)
    -a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z
end

#Exponentiation
function expb(a::MVeven)
    a = project(a,2)
    nrm = sqrt(dot(a,-a))
    if iszero(nrm)
        return MVeven(one(nrm),0,0,0)
    else
        return cos(nrm) + sin(nrm)*a/nrm
    end
end

function exp(a::MVeven)
    R = expb(a)
    if iszero(a.w)
        return R
    else 
        return exp(a.w)*R
    end
end

#Comparison. Using default tolerances.
Base.isapprox(a::MVeven, b::MVeven) = isapprox(a.w,b.w) && isapprox(a.x,b.x) && isapprox(a.y,b.y) && isapprox(a.z,b.z) 
Base.isapprox(a::MVodd, b::MVodd) = isapprox(a.w,b.w) && isapprox(a.x,b.x) && isapprox(a.y,b.y) && isapprox(a.z,b.z) 
