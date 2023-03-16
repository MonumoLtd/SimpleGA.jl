#=
Core code for the implementation of GA(4,1).
Representation is as a 2x2 matrix of quaternions.
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
    q1::Quaternion
    q2::Quaternion
    q3::Quaternion
    q4::Quaternion
end

struct MVodd
    q1::Quaternion
    q2::Quaternion
    q3::Quaternion
    q4::Quaternion
end

const qone = Quaternion(1,0,0,0)

#Addition / subtraction
function -(a::MVeven)
    MVeven(-a.q1,-a.q2,-a.q3,-a.q4)
end

function -(a::MVodd)
    MVodd(-a.q1,-a.q2,-a.q3,-a.q4)
end

function +(a::MVeven,b::MVeven)
    MVeven(a.q1 + b.q1, a.q2 + b.q2, a.q3 + b.q3, a.q4 + b.q4)
end

function +(a::MVodd,b::MVodd)
    MVodd(a.q1 + b.q1, a.q2 + b.q2, a.q3 + b.q3, a.q4 + b.q4)
end

function -(a::MVeven,b::MVeven)
    MVeven(a.q1 - b.q1, a.q2 - b.q2, a.q3 - b.q3, a.q4 - b.q4)
end

function -(a::MVodd,b::MVodd)
    MVodd(a.q1 - b.q1, a.q2 - b.q2, a.q3 - b.q3, a.q4 - b.q4)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::MVeven)
    MVeven(a.q1 + convert(Float64,num), a.q2, a.q3, a.q4 + convert(Float64,num))
end

function -(num::Number,a::MVeven)
    MVeven(-a.q1 + convert(Float64,num), -a.q2, -a.q3, -a.q4 + convert(Float64,num))
end

#Multiplication
function *(num::Number,a::MVeven)
    num = convert(Float64,num)
    MVeven(num*a.q1, num*a.q2, num*a.q3, num*a.q4)
end

function *(num::Number,a::MVodd)
    num = convert(Float64,num)
    MVodd(num*a.q1, num*a.q2, num*a.q3, num*a.q4)
end

function *(a::MVeven,b::MVeven)
    MVeven( a.q1*b.q1 - a.q2*b.q3,
            a.q1*b.q2 + a.q2*b.q4,
            a.q3*b.q1 + a.q4*b.q3,
            a.q4*b.q4 - a.q3*b.q2)    
end

function *(a::MVeven,b::MVodd)
    MVodd( a.q1*b.q1 - a.q2*b.q3,
            a.q1*b.q2 + a.q2*b.q4,
            a.q3*b.q1 + a.q4*b.q3,
            a.q4*b.q4 - a.q3*b.q2)    
end

function *(a::MVodd,b::MVeven)
    MVodd( a.q1*b.q1 - a.q2*b.q3,
            a.q1*b.q2 + a.q2*b.q4,
            a.q3*b.q1 + a.q4*b.q3,
            a.q4*b.q4 - a.q3*b.q2)    
end

function *(a::MVodd,b::MVodd)
    MVeven( -a.q1*b.q1 + a.q2*b.q3,
            -a.q1*b.q2 - a.q2*b.q4,
            -a.q3*b.q1 - a.q4*b.q3,
            -a.q4*b.q4 + a.q3*b.q2)    
end

#Division by a real
function /(a::MVeven,num::Number)
    num = convert(Float64,1/num)
    MVeven(num*a.q1, num*a.q2, num*a.q3, num*a.q4)
end

function /(a::MVodd,num::Number)
    num = convert(Float64,1/num)
    MVodd(num*a.q1, num*a.q2, num*a.q3, num*a.q4)
end

#Reverse
function adjoint(a::MVeven)
    MVeven(conj(a.q4),conj(a.q2), conj(a.q3), conj(a.q1))
end

function adjoint(a::MVodd)
    MVodd(conj(a.q4),conj(a.q2), conj(a.q3), conj(a.q1))
end

#Grade and projection
function project(a::MVeven,n::Integer)
    if (n==0)
        return MVeven(0.5*(a.q1.w+a.q4.w)*qone, qzero, qzero, 0.5*(a.q1.w+a.q4.w)*qone)
    elseif (n==2)
        qtmp = imag_part(0.5*(a.q1+a.q4))
        stmp = 0.5*(a.q1.w - a.q4.w)
        return MVeven(stmp+qtmp, imag_part(a.q2), imag_part(a.q3), -stmp+ qtmp )
    elseif (n==4)
        qtmp = imag_part(0.5*(a.q1-a.q4))
        return MVeven(qtmp, real_part(a.q2), real_part(a.q3), -qtmp )
    else
        return MVeven(qzero,qzero,qzero,qzero)
    end
end

function project(a::MVodd,n::Integer)
    if (n==5)
        return MVodd(0.5*(a.q1.w+a.q4.w)*qone, qzero, qzero, 0.5*(a.q1.w+a.q4.w)*qone)
    elseif (n==3)
        qtmp = imag_part(0.5*(a.q1+a.q4))
        stmp = 0.5*(a.q1.w - a.q4.w)
        return MVodd(stmp+qtmp, imag_part(a.q2), imag_part(a.q3), -stmp+ qtmp )
    elseif (n==1)
        qtmp = imag_part(0.5*(a.q1-a.q4))
        return MVodd(qtmp, real_part(a.q2), real_part(a.q3), -qtmp )
    else
        return MVodd(qzero,qzero,qzero,qzero)
    end
end

function tr(a::MVeven)
    0.5*(a.q1.w + a.q4.w)
end

function dot(a::MVeven, b::MVeven)
    0.5*(dot(a.q1,b.q1) - dot(a.q2,b.q3) + dot(a.q4,b.q4) - dot(a.q3,b.q2))    
end

function dot(a::MVodd, b::MVodd)
    -0.5*(dot(a.q1,b.q1) - dot(a.q2,b.q3) + dot(a.q4,b.q4) - dot(a.q3,b.q2))    
end

#Exponentiation
#Uses a simple scale and square implementation.
function exp(a::MVeven)
    anrm = dot(a.q1,conj(a.q1))+dot(a.q2,conj(a.q2))+dot(a.q3,conj(a.q3))+dot(a.q4,conj(a.q4))
    s = max(ceil(Int,log(2,anrm))-1,0)
    a = 1/2^s*a
    res = 1+a
    powa = a
    for i in 2:12
        powa *= a/i
        res += powa
    end
    while s > 0
        res = res*res
        s -= 1
    end
    return res
end

#Remove non-bivector terms. 
#TODO - investigate if closed form gives any performance benefits. Suspect all the if statements would slow this down.
function expb(a::MVeven)
    a = project(a,2)
    R = exp(a)
    delt = R*R'-1
    return (1-0.5*delt + 0.375*delt*delt)*R
end  


#Comparison
Base.isapprox(a::MVeven, b::MVeven) = isapprox(a.q1,b.q1) && isapprox(a.q2, b.q2) && isapprox(a.q3,b.q3) && isapprox(a.q4, b.q4)
Base.isapprox(a::MVodd, b::MVodd) = isapprox(a.q1,b.q1) && isapprox(a.q2, b.q2) && isapprox(a.q3,b.q3) && isapprox(a.q4, b.q4) 

