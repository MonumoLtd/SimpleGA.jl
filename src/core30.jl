 #=
Core code for the implementation of GA(3,0).
Underlying representation is with quaternions, though in this case we do not use the quaternion code as that would put in an unnecesary layer of indirection.
Instead we just use the hard-coded version of quaternions multiplication.
=#


import ..project
import ..expb

struct Even{T<:Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

struct Odd{T<:Real} <: Number
    w::T
    x::T
    y::T
    z::T
end

Base.convert(::Type{Even{T}},a::Even) where {T <: Real} = Even{T}(convert(T,a.w), convert(T,a.x), convert(T,a.y), convert(T,a.z))
Base.convert(::Type{Odd{T}},a::Odd) where {T <: Real} = Odd{T}(convert(T,a.w), convert(T,a.x), convert(T,a.y), convert(T,a.z))


#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.w,-a.x,-a.y,-a.z)
Base.:(-)(a::Odd) = Odd(-a.w,-a.x,-a.y,-a.z)
Base.:(+)(a::Even,b::Even) = Even(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
Base.:(+)(a::Odd,b::Odd) = Odd(a.w + b.w, a.x + b.x, a.y + b.y, a.z + b.z)
Base.:(-)(a::Even,b::Even) = Even(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)
Base.:(-)(a::Odd,b::Odd) = Odd(a.w - b.w, a.x - b.x, a.y - b.y, a.z - b.z)


#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number,a::Even) = Even(a.w+num, a.x, a.y, a.z)
Base.:(-)(num::Number,a::Even) = Even(-a.w+num, -a.x, -a.y, -a.z)
   

#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.w, num*a.x, num*a.y, num*a.z)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.w, num*a.x, num*a.y, num*a.z)

function Base.:(*)(a::Even, b::Even)
    Even(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
             a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
             a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
             a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function Base.:(*)(a::Even, b::Odd)
    Odd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function Base.:(*)(a::Odd, b::Even)
    Odd(a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y + a.y*b.w + a.z*b.x - a.x*b.z,
            a.w*b.z + a.z*b.w + a.x*b.y - a.y*b.x )
end

function Base.:(*)(a::Odd, b::Odd)
    Even(-a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z,
            -a.w*b.x - a.x*b.w - a.y*b.z + a.z*b.y,
            -a.w*b.y - a.y*b.w - a.z*b.x + a.x*b.z,
            -a.w*b.z - a.z*b.w - a.x*b.y + a.y*b.x )
end


#Reverse
LinearAlgebra.adjoint(a::Even) = Even(a.w,-a.x, -a.y, -a.z)
LinearAlgebra.adjoint(a::Odd) = Odd(-a.w,a.x, a.y, a.z)


#Grade and projection
function project(a::Even,n::Integer)
    if (n == 0)
        return Even(a.w,zero(a.w),zero(a.w),zero(a.w))
    elseif (n == 2)
        return Even(zero(a.w),a.x,a.y,a.z)
    else
        return Even(zero(a.w),zero(a.w),zero(a.w),zero(a.w))
    end
end


function project(a::Odd,n::Integer)
    if (n == 3)
        return Odd(a.w,zero(a.w),zero(a.w),zero(a.w))
    elseif (n == 1)
        return Odd(zero(a.w),a.x,a.y,a.z)
    else
        return Odd(zero(a.w),zero(a.w),zero(a.w),zero(a.w))
    end
end


LinearAlgebra.tr(a::Even) = a.w
LinearAlgebra.dot(a::Even, b::Even) = a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
LinearAlgebra.dot(a::Odd, b::Odd) = -a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z


#Exponentiation
function expb(a::Even)
    a = project(a,2)
    nrm = sqrt(dot(a,-a))
    if iszero(nrm)
        return Even(one(a.w),zero(a.w),zero(a.w),zero(a.w))
    else
        return cos(nrm) + sin(nrm)*a/nrm
    end
end

function Base.exp(a::Even)
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
