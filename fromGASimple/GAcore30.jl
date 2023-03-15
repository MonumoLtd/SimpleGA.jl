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

struct Even
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

struct Odd
    w::Float64
    x::Float64
    y::Float64
    z::Float64
end

#Addition / subtraction
function -(a::Even)
    Even(-a.w,-a.x,-a.y,-a.z)
end

function -(a::Odd)
    Odd(-a.w,-a.x,-a.y,-a.z)
end

function +(a::Even, b::Even)
    Even(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function (-)(a::Even,b::Even)
    Even(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

function +(a::Odd, b::Odd)
    Odd(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
end

function -(a::Odd, b::Odd)
    Odd(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::Even)
    Even(a.w+convert(Float64,num), a.x, a.y, a.z)
end

function -(num::Number,a::Even)
    Even(-a.w+convert(Float64,num), -a.x, -a.y, -a.z)
end

#Multiplication
function *(num::Number,a::Even)
    num = convert(Float64,num)
    Even(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(num::Number,a::Odd)
    num = convert(Float64,num)
    Odd(num*a.w, num*a.x, num*a.y, num*a.z)
end

function *(a::Even, b::Even)
    Even(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
             a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
             a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
             a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::Even, b::Odd)
    Odd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::Odd, b::Even)
    Odd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function *(a::Odd, b::Odd)
    Even(-a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z,
            -a.w*b.x - a.x*b.w - a.y*b.z + a.z*b.y,
            -a.w*b.y - a.y*b.w - a.z*b.x + a.x*b.z,
            -a.w*b.z - a.z*b.w - a.x*b.y + a.y*b.x )
end

#Division by a real
function /(a::Even,num::Number)
    num = convert(Float64,1/num)
    Even(num*a.w, num*a.x, num*a.y, num*a.z)
end

function /(a::Odd,num::Number)
    num = convert(Float64,1/num)
    Odd(num*a.w, num*a.x, num*a.y, num*a.z)
end

#Reverse
function adjoint(a::Even)
    Even(a.w,-a.x, -a.y, -a.z)
end

function adjoint(a::Odd)
    Odd(-a.w,a.x, a.y, a.z)
end

#Grade and projection
function project(a::Even,n::Integer)
    if (n == 0)
        return Even(a.w,0,0,0)
    elseif (n == 2)
        return Even(0,a.x,a.y,a.z)
    else
        return Even(0,0,0,0)
    end
end

function project(a::Odd,n::Integer)
    if (n == 3)
        return Odd(a.w,0,0,0)
    elseif (n == 1)
        return Odd(0,a.x,a.y,a.z)
    else
        return Odd(0,0,0,0)
    end
end

function tr(a::Even)
    a.w
end

function dot(a::Even, b::Even)
    a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
end

function dot(a::Odd, b::Odd)
    -a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z
end

#Exponentiation
function expb(a::Even)
    a = project(a,2)
    nrm = sqrt(dot(a,-a))
    if iszero(nrm)
        return Even(one(nrm),0,0,0)
    else
        return cos(nrm) + sin(nrm)*a/nrm
    end
end

function exp(a::Even)
    R = expb(a)
    if iszero(a.w)
        return R
    else 
        return exp(a.w)*R
    end
end

#Comparison. Using default tolerances.
Base.isapprox(a::Even, b::Even) = isapprox(a.w,b.w) && isapprox(a.x,b.x) && isapprox(a.y,b.y) && isapprox(a.z,b.z) 
Base.isapprox(a::Odd, b::Odd) = isapprox(a.w,b.w) && isapprox(a.x,b.x) && isapprox(a.y,b.y) && isapprox(a.z,b.z) 
