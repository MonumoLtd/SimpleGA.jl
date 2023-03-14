#=
Core code for the implementation of GA(4,0).
Work using self-dual and anti-self-dual decomposition, so just use a pair of quaternions.
Useful algebra for projective geometry.
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
    qp::Quaternion
    qm::Quaternion
end

struct MVodd
    qp::Quaternion
    qm::Quaternion
end

#Addition / subtraction
function -(a::MVeven)
    MVeven(-a.qp,-a.qm)
end

function -(a::MVodd)
    MVodd(-a.qp,-a.qm)
end

function +(a::MVeven,b::MVeven)
    MVeven(a.qp + b.qp, a.qm + b.qm)
end

function +(a::MVodd,b::MVodd)
    MVodd(a.qp + b.qp, a.qm + b.qm)
end

function -(a::MVeven,b::MVeven)
    MVeven(a.qp - b.qp, a.qm - b.qm)
end

function -(a::MVodd,b::MVodd)
    MVodd(a.qp - b.qp, a.qm - b.qm)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::MVeven)
    MVeven(a.qp + convert(Float64,num), a.qm + convert(Float64,num))
end

function -(num::Number,a::MVeven)
    MVeven(-a.qp + convert(Float64,num), -a.qm + convert(Float64,num))
end

#Multiplication
function *(num::Number,a::MVeven)
    num = convert(Float64,num)
    MVeven(num*a.qp,  num*a.qm)
end

function *(num::Number,a::MVodd)
    num = convert(Float64,num)
    MVodd(num*a.qp,  num*a.qm)
end

function *(a::MVeven,b::MVeven)
    MVeven( a.qp*b.qp,   a.qm*b.qm)  
end

function *(a::MVeven,b::MVodd)
    MVodd( a.qp*b.qp,   a.qm*b.qm)  
end

function *(a::MVodd,b::MVeven)
    MVodd( a.qp*b.qm, a.qm*b.qp)    
end

function *(a::MVodd,b::MVodd)
    MVeven( a.qp*b.qm, a.qm*b.qp)    
end

#Division by a real
function /(a::MVeven,num::Number)
    num = convert(Float64,1/num)
    MVeven(num*a.qp,  num*a.qm)
end

function /(a::MVodd,num::Number)
    num = convert(Float64,1/num)
    MVodd(num*a.qp,  num*a.qm)
end

#Reverse
function adjoint(a::MVeven)
    MVeven(conj(a.qp),conj(a.qm))
end

function adjoint(a::MVodd)
    MVodd(conj(a.qm),conj(a.qp))
end

#Grade and projection
function project(a::MVeven,n::Integer)
    if (n==0)
        return MVeven(Quaternion(0.5*(a.qp.w+a.qm.w),0,0,0), Quaternion(0.5*(a.qp.w+a.qm.w),0,0,0))
    elseif (n==2)
        return MVeven(imag_part(a.qp), imag_part(a.qm))
    elseif (n==4)
        return MVeven(Quaternion(0.5*(a.qp.w-a.qm.w),0,0,0), Quaternion(0.5*(-a.qp.w+a.qm.w),0,0,0))
    else
        return MVeven(qzero,qzero)
    end
end

function project(a::MVodd,n::Integer)
    if (n==1)
        return MVodd(0.5*(a.qp + conj(a.qm)), 0.5*(a.qm + conj(a.qp)))
    elseif (n==3)
        return MVodd(0.5*(a.qp - conj(a.qm)), 0.5*(a.qm - conj(a.qp)))
    else
        return MVodd(qzero, qzero)
    end
end

function tr(a::MVeven)
    0.5*(a.qp.w + a.qm.w)
end

function dot(a::MVeven, b::MVeven)
    0.5*(dot(a.qp,b.qp) + dot(a.qm,b.qm))
end

function dot(a::MVodd, b::MVodd)
    0.5*(dot(a.qp,b.qm) + dot(a.qm,b.qp))
end

#Exponentiation
function expb(a::MVeven)
    return MVeven(expb(a.qp), expb(a.qm))
end

function exp(a::MVeven)
    return MVeven(exp(a.qp), exp(a.qm))
end

#Comparison
Base.isapprox(a::MVeven, b::MVeven) = isapprox(a.qp,b.qp) && isapprox(a.qm, b.qm) 
Base.isapprox(a::MVodd, b::MVodd) = isapprox(a.qp,b.qp) && isapprox(a.qm, b.qm) 
