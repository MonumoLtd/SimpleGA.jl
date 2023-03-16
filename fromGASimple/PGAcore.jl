#=
Core code for the implementation of GA(3,0,1).
Work using a pair of quaternions. q is the main term and n is the null term.
Even / odd map performed by I3
Useful for performance if CGA is too slow.
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

using ..Quaternions

struct MVeven
    q::Quaternion
    n::Quaternion
end

struct MVodd
    q::Quaternion
    n::Quaternion
end

#Addition / subtraction
function -(a::MVeven)
    MVeven(-a.q,-a.n)
end

function -(a::MVodd)
    MVodd(-a.q,-a.n)
end

function +(a::MVeven,b::MVeven)
    MVeven(a.q + b.q, a.n + b.n)
end

function +(a::MVodd,b::MVodd)
    MVodd(a.q + b.q, a.n + b.n)
end

function -(a::MVeven,b::MVeven)
    MVeven(a.q - b.q, a.n - b.n)
end

function -(a::MVodd,b::MVodd)
    MVodd(a.q - b.q, a.n - b.n)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::MVeven)
    MVeven(a.q + convert(Float64,num), a.n)
end

function -(num::Number,a::MVeven)
    MVeven(-a.q + convert(Float64,num), -a.n )
end

#Multiplication
function *(num::Number,a::MVeven)
    num = convert(Float64,num)
    MVeven(num*a.q,  num*a.n)
end

function *(num::Number,a::MVodd)
    num = convert(Float64,num)
    MVodd(num*a.q,  num*a.n)
end

function *(a::MVeven,b::MVeven)
    MVeven( a.q*b.q, a.q*b.n + a.n*b.q)  
end

function *(a::MVeven,b::MVodd)
    MVodd( a.q*b.q, a.q*b.n - a.n*b.q)  
end

function *(a::MVodd,b::MVeven)
    MVodd( a.q*b.q, a.q*b.n + a.n*b.q)    
end

function *(a::MVodd,b::MVodd)
    MVeven( -a.q*b.q, -a.q*b.n + a.n*b.q )     
end

#Division by a real
function /(a::MVeven,num::Number)
    num = convert(Float64,1/num)
    MVeven(num*a.q,  num*a.n)
end

function /(a::MVodd,num::Number)
    num = convert(Float64,1/num)
    MVodd(num*a.q,  num*a.n)
end

#Reverse
function adjoint(a::MVeven)
    MVeven(conj(a.q),conj(a.n))
end

function adjoint(a::MVodd)
    MVodd(-conj(a.q),conj(a.n))
end

#Grade and projection
function project(a::MVeven,n::Integer)
    if (n==0)
        return MVeven(real_part(a.q), qzero)
    elseif (n==2)
        return MVeven(imag_part(a.q), imag_part(a.n))
    elseif (n==4)
        return MVeven(qzero, real_part(a.n))
    else
        return MVeven(qzero,qzero)
    end
end

function project(a::MVodd,n::Integer)
    if (n==1)
        return MVodd(imag_part(a.q), real_part(a.n))
    elseif (n==3)
        return MVodd(real_part(a.q), imag_part(a.n))
    else
        return MVodd(qzero, qzero)
    end
end

function tr(a::MVeven)
    a.q.w
end

function dot(a::MVeven, b::MVeven)
    dot(a.q,b.q)
end

function dot(a::MVodd, b::MVodd)
    -dot(a.q,b.q)
end

#Exponentiation
function expb(a::MVeven)
    a = project(a,2)
    aa = -a*a
    if iszero(tr(aa))
        return 1+a
    else
        f0 = sqrt(tr(aa))
        f1 = aa.n.w/2/f0
        cf = MVeven(Quaternion(cos(f0),0,0,0), Quaternion(-f1*sin(f0),0,0,0))
        sncf = MVeven(Quaternion(sin(f0)/f0,0,0,0),  Quaternion(f1/f0^2*(f0*cos(f0) - sin(f0)),0,0,0))
        return cf + sncf*a
    end
end

function exp(a::MVeven)
    R = expb(a)
    return exp(a.q.w)*(1+project(a,4))*R
end


#Comparison. Uses default tolerances.
Base.isapprox(a::MVeven, b::MVeven) = isapprox(a.q,b.q) && isapprox(a.n, b.n) 
Base.isapprox(a::MVodd, b::MVodd) = isapprox(a.q,b.q) && isapprox(a.n, b.n) 
