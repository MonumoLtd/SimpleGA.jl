#=
Core code for the implementation of GA(4,0).
Work using self-dual and anti-self-dual decomposition, so just use a pair of quaternions.
Useful algebra for projective geometry.
=#


import ..project
import ..expb

using ..Quaternions

struct Even{T<:Real} <: Number
    qp::Quaternion{T}
    qm::Quaternion{T}
end

struct Odd{T<:Real} <: Number
    qp::Quaternion{T}
    qm::Quaternion{T}
end

function Base.convert(::Type{Even{T}},a::Even) where {T <: Real} 
    return Even{T}(convert(Quaternion{T},a.qp), convert(Quaternion{T}, a.qm))
end

function Base.convert(::Type{Odd{T}},a::Odd) where {T <: Real} 
    return Odd{T}(convert(Quaternion{T},a.qp), convert(Quaternion{T}, a.qm))
end

Base.zero(a::Even) = Even(zero(a.qp), zero(a.qm))
Base.zero(a::Odd) = Odd(zero(a.qp), zero(a.qm))
Base.one(a::Even) = Even(one(a.qp), one(a.qm))

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.qp,-a.qm)
Base.:(-)(a::Odd) = Odd(-a.qp,-a.qm)
Base.:(+)(a::Even,b::Even) = Even(a.qp + b.qp, a.qm + b.qm)
Base.:(+)(a::Odd,b::Odd) = Odd(a.qp + b.qp, a.qm + b.qm)
Base.:(-)(a::Even,b::Even) = Even(a.qp - b.qp, a.qm - b.qm)
Base.:(-)(a::Odd,b::Odd) = Odd(a.qp - b.qp, a.qm - b.qm)


#Scalar addition / subtraction. Other cases are in GAcommon
Base.:(+)(num::Number,a::Even) = Even(a.qp + num, a.qm + num)
Base.:(-)(num::Number,a::Even) = Even(-a.qp +num, -a.qm + num)

#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.qp,  num*a.qm)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.qp,  num*a.qm)
Base.:(*)(a::Even,b::Even) = Even( a.qp*b.qp,   a.qm*b.qm)
Base.:(*)(a::Even,b::Odd) = Odd( a.qp*b.qp,   a.qm*b.qm) 
Base.:(*)(a::Odd, b::Even) = Odd( a.qp*b.qm, a.qm*b.qp)     
Base.:(*)(a::Odd, b::Odd) = Even( a.qp*b.qm, a.qm*b.qp)     


#Reverse
LinearAlgebra.adjoint(a::Even) = Even(conj(a.qp),conj(a.qm))
LinearAlgebra.adjoint(a::Odd) = Odd(conj(a.qm),conj(a.qp))


#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        return Even(Quaternion((a.qp.w+a.qm.w)/2), Quaternion((a.qp.w+a.qm.w)/2))
    elseif (n==2)
        return (a-a')/2
    elseif (n==4)
        return Even(Quaternion((a.qp.w-a.qm.w)/2), Quaternion((-a.qp.w+a.qm.w)/2))
    else
        return zero(a)
    end
end

function project(a::Odd,n::Integer)
    if (n==1)
        return Odd((a.qp + conj(a.qm))/2, (a.qm + conj(a.qp))/2)
    elseif (n==3)
        return Odd((a.qp - conj(a.qm))/2, (a.qm - conj(a.qp))/2)
    else
        return zero(a)
    end
end

LinearAlgebra.tr(a::Even) = (a.qp.w + a.qm.w)/2
LinearAlgebra.dot(a::Even, b::Even) = (dot(a.qp,b.qp) + dot(a.qm,b.qm))/2
LinearAlgebra.dot(a::Odd, b::Odd) = (dot(a.qp,b.qm) + dot(a.qm,b.qp))/2


#Exponentiation
function expb(a::Even)
    return Even(expb(a.qp), expb(a.qm))
end

function Base.exp(a::Even)
    return Even(exp(a.qp), exp(a.qm))
end

#Comparison
Base.isapprox(a::Even, b::Even) = isapprox(a.qp,b.qp) && isapprox(a.qm, b.qm) 
Base.isapprox(a::Odd, b::Odd) = isapprox(a.qp,b.qp) && isapprox(a.qm, b.qm) 
