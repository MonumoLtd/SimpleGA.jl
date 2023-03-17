#=
Core code for the implementation of GA(3,0,1).
Work using a pair of quaternions. q is the main term and n is the null term.
Even / odd map performed by I3
Useful for performance if CGA is too slow.
=#


import ..project
import ..expb

using ..Quaternions

struct Even{T<:Real} <: Number
    q::Quaternion{T}
    n::Quaternion{T}
end

struct Odd{T<:Real} <: Number
    q::Quaternion{T}
    n::Quaternion{T}
end


function Base.convert(::Type{Even{T}},a::Even) where {T <: Real} 
    return Even{T}(convert(Quaternion{T},a.q), convert(Quaternion{T}, a.n))
end

function Base.convert(::Type{Odd{T}},a::Odd) where {T <: Real} 
    return Odd{T}(convert(Quaternion{T},a.q), convert(Quaternion{T}, a.n))
end

Base.zero(a::Even) = Even(zero(a.q), zero(a.n))
Base.zero(a::Odd) = Odd(zero(a.q), zero(a.n))
Base.one(a::Even) = Even(one(a.q), zero(a.n))


#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.q,-a.n)
Base.:(-)(a::Odd) = Odd(-a.q,-a.n)
Base.:(+)(a::Even,b::Even) = Even(a.q + b.q, a.n + b.n)
Base.:(+)(a::Odd,b::Odd) = Odd(a.q + b.q, a.n + b.n)
Base.:(-)(a::Even,b::Even) = Even(a.q - b.q, a.n - b.n)
Base.:(-)(a::Odd,b::Odd) = Odd(a.q - b.q, a.n - b.n)

#Scalar addition / subtraction. Other cases are in GAcommon
Base.:(+)(num::Number,a::Even) = Even(a.q + num, a.n)
Base.:(-)(num::Number,a::Even) = Even(-a.q +num, -a.qm + num)

#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.q,  num*a.n)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.q,  num*a.n)
Base.:(*)(a::Even,b::Even) = Even( a.q*b.q, a.q*b.n + a.n*b.q)
Base.:(*)(a::Even,b::Odd) = Odd( a.q*b.q, a.q*b.n - a.n*b.q) 
Base.:(*)(a::Odd, b::Even) = Odd( a.q*b.q, a.q*b.n + a.n*b.q)     
Base.:(*)(a::Odd, b::Odd) = Even( -a.q*b.q, -a.q*b.n + a.n*b.q)    


#Reverse
LinearAlgebra.adjoint(a::Even) = Even(conj(a.q),conj(a.n))
LinearAlgebra.adjoint(a::Odd) = Odd(-conj(a.q),conj(a.n))


#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        return Even(real_part(a.q), zero(a.q))
    elseif (n==2)
        return Even(imag_part(a.q), imag_part(a.n))
    elseif (n==4)
        return Even(zero(a.q), real_part(a.n))
    else
        return zero(a)
    end
end

function project(a::Odd,n::Integer)
    if (n==1)
        return Odd(imag_part(a.q), real_part(a.n))
    elseif (n==3)
        return Odd(real_part(a.q), imag_part(a.n))
    else
        return zero(a)
    end
end


LinearAlgebra.tr(a::Even) = a.q.w 
LinearAlgebra.dot(a::Even, b::Even) = dot(a.q,b.q)
LinearAlgebra.dot(a::Odd, b::Odd) = -dot(a.q,b.q)


#Exponentiation
function expb(a::Even)
    a = project(a,2)
    aa = -a*a
    if iszero(tr(aa))
        return 1+a
    else
        f0 = sqrt(tr(aa))
        f1 = aa.n.w/2/f0
        cf = Even(Quaternion(cos(f0),0,0,0), Quaternion(-f1*sin(f0),0,0,0))
        sncf = Even(Quaternion(sin(f0)/f0,0,0,0),  Quaternion(f1/f0^2*(f0*cos(f0) - sin(f0)),0,0,0))
        return cf + sncf*a
    end
end

function Base.exp(a::Even)
    R = expb(a)
    return exp(a.q.w)*(1+project(a,4))*R
end


#Comparison. Uses default tolerances.
Base.isapprox(a::Even, b::Even) = isapprox(a.q,b.q) && isapprox(a.n, b.n) 
Base.isapprox(a::Odd, b::Odd) = isapprox(a.q,b.q) && isapprox(a.n, b.n) 
Base.isequal(a::Even, b::Even) = isequal(a.q,b.q) && isequal(a.n, b.n) 
Base.isequal(a::Odd, b::Odd) = isequal(a.q,b.q) && isequal(a.n, b.n) 
