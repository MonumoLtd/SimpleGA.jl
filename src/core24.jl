#=
Core code for the implementation of GA(2,4).
Base element is a 4x4 Complex matrix built on Static Arrays library.
This is the conformal algebra for spacetime, also relevant to twistor geometry.
=#

import ..project
import ..expb

struct Even{T<:Real} <: Number
    m::SMatrix{4,4,Complex{T},16}
end

struct Odd{T<:Real} <: Number
    m::SMatrix{4,4,Complex{T},16}
end

const id4 = SMatrix{4,4,Complex{Int8},16}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
const adj = SMatrix{4,4,Complex{Int8},16}([0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ])
const rev = SMatrix{4,4,Complex{Int8},16}([0  0 0 -im; 0 0 im 0;  0 -im 0 0; im 0 0 0])



function Base.convert(::Type{Even{T}},a::Even) where {T <: Real} 
    return Even{T}(convert(SMatrix{4,4,Complex{T},16}, (a.m) ))
end

function Base.convert(::Type{Odd{T}},a::Even) where {T <: Real} 
    return Odd{T}(convert(SMatrix{4,4,Complex{T},16}, a.m) )
end

Base.zero(a::Even) = Even(zero(a.m))
Base.zero(a::Odd) = Odd(zero(a.m))
Base.one(a::Even) = Even(one(a.m))


#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.m)
Base.:(-)(a::Odd) = Odd(-a.m)
Base.:(+)(a::Even,b::Even) = Even(a.m + b.m)
Base.:(+)(a::Odd,b::Odd) = Odd(a.m + b.m)
Base.:(-)(a::Even,b::Even) =  Even(a.m - b.m)
Base.:(-)(a::Odd,b::Odd) = Odd(a.m - b.m)


#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number,a::Even) = Even(a.m + num*one(a.m))
Base.:(-)(num::Number,a::Even) = Even(-a.m + num*one(a.m))


#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.m)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.m)
Base.:(*)(a::Even,b::Even) = Even(a.m*b.m)  
Base.:(*)(a::Even,b::Odd) = Odd(a.m*b.m)  
Base.:(*)(a::Odd,b::Even) = Odd( a.m*g2.m*conj(b.m)*g2.m)    
Base.:(*)(a::Odd,b::Odd) = Even( a.m*g2.m*conj(b.m)*g2.m)    


#Reverse
LinearAlgebra.adjoint(a::Even) = Even(-adj*(a.m)'*adj)
LinearAlgebra.adjoint(a::Odd) = Odd(rev*transpose(a.m)*rev)


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
        return zero(a)
    end
end



function project(a::Odd,n::Integer)
    if (n==3)
        return (a - a')/2
    elseif (n==1)
        tmp = (a+a')*g0/2
        return (project(tmp,0)+project(tmp,2))*g0
    elseif (n==5)
        tmp = (a+a')*g0/2
        return (project(tmp,4)+project(tmp,6))*g0
    else
        return zero(a)
    end
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
