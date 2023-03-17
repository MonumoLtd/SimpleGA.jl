#=
Core code for the implementation of GA(3,3).
Work using self-dual and anti-self-dual decomposition. Base element is a 4x4 matrix built on Static Arrays library.
Useful algebra for line geometry.
=#

import ..project
import ..expb


struct Even{T<:Real} <: Number
    p::SMatrix{4,4,T,16}
    m::SMatrix{4,4,T,16}
end

struct Odd{T<:Real} <: Number
    p::SMatrix{4,4,T,16}
    m::SMatrix{4,4,T,16}
end

id4 = SMatrix{4,4,Int8}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
J4 = SMatrix{4,4,Int8}([0 -1 0 0; 1 0 0 0; 0 0 0 1; 0 0 -1 0])

function Base.convert(::Type{Even{T}},a::Even) where {T <: Real} 
    return Even{T}(convert(SMatrix{4,4,T,16}, a.p),convert(SMatrix{4,4,T,16}, a.m) )
end

function Base.convert(::Type{Odd{T}},a::Even) where {T <: Real} 
    return Odd{T}(convert(SMatrix{4,4,T,16}, a.p),convert(SMatrix{4,4,T,16}, a.m) )
end

Base.zero(a::Even) = Even(zero(a.p), zero(a.m))
Base.zero(a::Odd) = Odd(zero(a.p), zero(a.m))
Base.one(a::Even) = Even(one(a.p), one(a.m))


#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.p,-a.m)
Base.:(-)(a::Odd) = Odd(-a.p,-a.m)
Base.:(+)(a::Even,b::Even) = Even(a.p + b.p, a.m + b.m)
Base.:(+)(a::Odd,b::Odd) = Odd(a.p + b.p, a.m + b.m)
Base.:(-)(a::Even,b::Even) =  Even(a.p - b.p, a.m - b.m)
Base.:(-)(a::Odd,b::Odd) = Odd(a.p - b.p, a.m - b.m)



#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number,a::Even) = Even(a.p + num*one(a.p),  a.m + num*one(a.m))
Base.:(-)(num::Number,a::Even) = Even(-a.p + num*one(a.p),  -a.m + num*one(a.m))


#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.p,  num*a.m)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.p,  num*a.m)
Base.:(*)(a::Even,b::Even) = Even( a.p*b.p,   a.m*b.m)  
Base.:(*)(a::Even,b::Odd) = Odd( a.p*b.p,   a.m*b.m)  
Base.:(*)(a::Odd,b::Even) = Odd( a.p*b.m, a.m*b.p)    
Base.:(*)(a::Odd,b::Odd) = Even( a.p*b.m, a.m*b.p)    


#Reverse
function LinearAlgebra.adjoint(a::Even) 
    Even(-J4*transpose(a.m)*J4, -J4*transpose(a.p)*J4 )
end

function LinearAlgebra.adjoint(a::Odd)
    Odd(-J4*transpose(a.p)*J4, -J4*transpose(a.m)*J4 )
end


#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        scl = (tr(a.p)+tr(a.m))/8
        return scl*one(a)
    elseif (n==2)
        scl = (tr(a.p)-tr(a.m))/8
        return (a-a')/2 - scl*Even(one(a.p),-one(a.m))
    elseif (n==4)
        scl = (tr(a.p)+tr(a.m))/8
        return (a+a')/2 - scl*one(a)
    elseif (n == 6)
        scl = (tr(a.p)-tr(a.m))/8
        return scl*Even(one(a.p),-one(a.m))
    else
        return zero(a)
    end
end

function project(a::Odd,n::Integer)
    if (n==3)
        return (a - a')/2
    elseif (n==1)
        tmp = (a+a')*e3/2
        return (project(tmp,0)+project(tmp,2))*e3
    elseif (n==5)
        tmp = (a+a')*e3/2
        return (project(tmp,4)+project(tmp,6))*e3
    else
        return zero(a)
    end
end   


LinearAlgebra.tr(a::Even) = (tr(a.p)+tr(a.m))/8
LinearAlgebra.dot(a::Even, b::Even) = (tr(a.p*b.p)+tr(a.m*b.m))/8   
LinearAlgebra.dot(a::Odd, b::Odd) = (tr(a.p*b.m)+tr(a.m*b.p))/8   


#Exponentiation
function Base.exp(a::Even)
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
