#=
Core code for the implementation of GA(2,4).
Base element is a 4x4 complex matrix built on Static Arrays library.
This is the conformal algebra for spacetime, also relevant to twistor geometry.
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
    m::SArray{Tuple{4,4},ComplexF64,2,16}
end

struct Odd
    m::SArray{Tuple{4,4},ComplexF64,2,16}
end

complexSA(arr) = SMatrix{4,4,ComplexF64,16}(convert(Matrix{ComplexF64},arr))

const id4 = complexSA([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
const adj = complexSA([0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ])
const rev = complexSA([0  0 0 -im; 0 0 im 0;  0 -im 0 0; im 0 0 0])
const mzero = complexSA([0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0])
 

#Addition / subtraction
function -(a::Even)
    Even(-a.m)
end

function -(a::Odd)
    Odd(-a.m)
end

function +(a::Even,b::Even)
    Even(a.m + b.m)
end

function +(a::Odd,b::Odd)
    Odd(a.m + b.m)
end

function -(a::Even,b::Even)
    Even(a.m - b.m)
end

function -(a::Odd,b::Odd)
    Odd(a.m - b.m)
end


#Scalar addition / subtraction. Other cases are in GAcommon
function +(num::Number,a::Even)
    Even(a.m + convert(ComplexF64,num)*id4)
end

function -(num::Number,a::Even)
    Even(-a.m + convert(ComplexF64,num)*id4)
end

function *(num::Number,a::Even)
    num = convert(Float64,num)
    Even(num*a.m)
end

function *(num::Number,a::Odd)
    num = convert(Float64,num)
    Odd(num*a.m)
end

function *(a::Even,b::Even)
    Even(a.m*b.m)  
end

function *(a::Even,b::Odd)
    Odd( a.m*b.m)  
end

function *(a::Odd,b::Even)
    Odd( a.m*g2.m*conj(b.m)*g2.m)    
end

function *(a::Odd,b::Odd)
    Even(a.m*g2.m*conj(b.m)*g2.m)    
end


#Division by a real
function /(a::Even,num::Number)
    num = convert(Float64,1/num)
    Even(num*a.m)
end

function /(a::Odd,num::Number)
    num = convert(Float64,1/num)
    Odd(num*a.m)
end

#Reverse
function adjoint(a::Even)
    Even(-adj*(a.m)'*adj)
end

function adjoint(a::Odd)
    Odd(rev*transpose(a.m)*rev)
end


#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        scl = real((tr(a.m))/4)
        return Even(scl*id4)
    elseif (n==2)
        tmp = (a-a')/2
        return tmp - Even(tr(tmp.m)/4*id4)
    elseif (n==4)
        tmp = (a+a')/2
        return tmp - Even(tr(tmp.m)/4*id4)
    elseif (n == 6)
        scl = im*imag((tr(a.m))/4)
        return Even(scl*id4)
    else
        return Even(mzero)
    end
end


function project(a::Odd,n::Integer)
    if (n==3)
        return 0.5*(a - a')
    elseif (n==1)
        tmp = (a+a')*g0/2
        return (project(tmp,0)+project(tmp,2))*g0
    elseif (n==5)
        tmp = (a+a')*g0/2
        return (project(tmp,4)+project(tmp,6))*g0
    else
        return Odd(mzero)
    end
end   

function tr(a::Even)
    real(tr(a.m))/4
end

function dot(a::Even, b::Even)
    real(tr(a.m*b.m))/4  
end


function dot(a::Odd, b::Odd)
    bp = g2.m*conj(b.m)*g2.m
    real(tr(a.m*bp))/4  
end



#Exponentiation
function exp(a::Even)
    return Even(exp(a.m))
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
Base.isapprox(a::Even, b::Even) =  isapprox(a.m,b.m)
Base.isapprox(a::Odd, b::Odd) =  isapprox(a.m,b.m)
