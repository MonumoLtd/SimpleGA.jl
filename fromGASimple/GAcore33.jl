#=
Core code for the implementation of GA(3,3).
Work using self-dual and anti-self-dual decomposition. Base element is a 4x4 matrix built on Static Arrays library.
Useful algebra for line geometry.
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
    p::SArray{Tuple{4,4},Float64,2,16}
    m::SArray{Tuple{4,4},Float64,2,16}
end

struct Odd
    p::SArray{Tuple{4,4},Float64,2,16}
    m::SArray{Tuple{4,4},Float64,2,16}
end

const id4 = SA_F64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const J4 = SA_F64[0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0]
const mzero = SA_F64[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]


#Addition / subtraction
function -(a::Even)
    Even(-a.p,-a.m)
end

function -(a::Odd)
    Odd(-a.p,-a.m)
end

function +(a::Even,b::Even)
    Even(a.p + b.p, a.m + b.m)
end

function +(a::Odd,b::Odd)
    Odd(a.p + b.p, a.m + b.m)
end

function -(a::Even,b::Even)
    Even(a.p - b.p, a.m - b.m)
end

function -(a::Odd,b::Odd)
    Odd(a.p - b.p, a.m - b.m)
end

#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::Even)
    Even(a.p + convert(Float64,num)*id4, a.m + convert(Float64,num)*id4)
end

function -(num::Number,a::Even)
    Even(-a.p + convert(Float64,num)*id4, -a.m + convert(Float64,num)*id4)
end

function *(num::Number,a::Even)
    num = convert(Float64,num)
    Even(num*a.p,  num*a.m)
end

function *(num::Number,a::Odd)
    num = convert(Float64,num)
    Odd(num*a.p,  num*a.m)
end

function *(a::Even,b::Even)
    Even( a.p*b.p,   a.m*b.m)  
end

function *(a::Even,b::Odd)
    Odd( a.p*b.p,   a.m*b.m)  
end

function *(a::Odd,b::Even)
    Odd( a.p*b.m, a.m*b.p)    
end

function *(a::Odd,b::Odd)
    Even( a.p*b.m, a.m*b.p)    
end

#Division by a real
function /(a::Even,num::Number)
    num = convert(Float64,1/num)
    Even(num*a.p,  num*a.m)
end

function /(a::Odd,num::Number)
    num = convert(Float64,1/num)
    Odd(num*a.p,  num*a.m)
end

#Reverse
function adjoint(a::Even)
    Even(-J4*transpose(a.m)*J4, -J4*transpose(a.p)*J4 )
end

function adjoint(a::Odd)
    Odd(-J4*transpose(a.p)*J4, -J4*transpose(a.m)*J4 )
end

#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        scl = (tr(a.p)+tr(a.m))/8
        return Even(scl*id4,scl*id4)
    elseif (n==2)
        scl = (tr(a.p)-tr(a.m))/8
        return 0.5*(a - a') - Even(scl*id4,-scl*id4)
    elseif (n==4)
        scl = (tr(a.p)+tr(a.m))/8
        return 0.5*(a+a') - Even(scl*id4,scl*id4)
    elseif (n == 6)
        scl = (tr(a.p)-tr(a.m))/8
        return Even(scl*id4,-scl*id4)
    else
        return Even(mzero,mzero)
    end
end

function project(a::Odd,n::Integer)
    if (n==3)
        return 0.5*(a - a')
    elseif (n==1)
        tmp = (a+a')*e3/2
        return (project(tmp,0)+project(tmp,2))*e3
    elseif (n==5)
        tmp = (a+a')*e3/2
        return (project(tmp,4)+project(tmp,6))*e3
    else
        return Odd(mzero,mzero)
    end
end   

function tr(a::Even)
    (tr(a.p)+tr(a.m))/8
end

#StaticArrays seems to have optimised the trace of a product operation already.
function dot(a::Even, b::Even)
    (tr(a.p*b.p)+tr(a.m*b.m))/8   
end

function dot(a::Odd, b::Odd)
    (tr(a.p*b.m)+tr(a.m*b.p))/8   
end


#Exponentiation
function exp(a::Even)
    return Even(exp(a.p),exp(a.m))
end

#TODO Any improvement here?
function expb(a::Even)
    a = project(a,2)
    R = exp(a)
    delt = R*R'-1
    return (1-0.5*delt + 0.375*delt*delt)*R
end

#Comparison
#StaticArrays does seem to lose some accuracy.
Base.isapprox(a::Even, b::Even) = isapprox(a.p,b.p) && isapprox(a.m,b.m)
Base.isapprox(a::Odd, b::Odd) = isapprox(a.p,b.p) && isapprox(a.m,b.m)
