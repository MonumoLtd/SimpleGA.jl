#=
Core code for the implementation of GA(1,3).
Underlying representation is with 2x2 complex matrices, though no matrices are formally constructed. The multiplication rules are hard-coded for efficiency.
Makes use of Julia's internal ComplexF64 format.
=#


import ..project
import ..expb


struct Even{T<:Real} <: Number
    c1::Complex{T}
    c2::Complex{T}
    c3::Complex{T}
    c4::Complex{T}
end

struct Odd{T<:Real} <: Number
    c1::Complex{T}
    c2::Complex{T}
    c3::Complex{T}
    c4::Complex{T}
end

function Base.convert(::Type{Even{T}},a::Even) where {T <: Real} 
    return Even{T}(convert(Complex{T},a.c1), convert(Complex{T},a.c2), convert(Complex{T},a.c3), convert(Complex{T},a.c4))
end

function Base.convert(::Type{Odd{T}},a::Odd) where {T <: Real} 
    return Odd{T}(convert(Complex{T},a.c1), convert(Complex{T},a.c2), convert(Complex{T},a.c3), convert(Complex{T},a.c4))
end

#Addition / subtraction
Base.:(-)(a::Even) = Even(-a.c1,-a.c2,-a.c3,-a.c4)
Base.:(-)(a::Odd) = Odd(-a.c1,-a.c2,-a.c3,-a.c4)
Base.:(+)(a::Even,b::Even) = Even(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(+)(a::Odd,b::Odd) = Odd(a.c1 + b.c1, a.c2 + b.c2, a.c3 + b.c3, a.c4 + b.c4)
Base.:(-)(a::Even,b::Even) =  Even(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)
Base.:(-)(a::Odd,b::Odd) = Odd(a.c1 - b.c1, a.c2 - b.c2, a.c3 - b.c3, a.c4 - b.c4)


#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number,a::Even) = Even(a.c1 + num, a.c2, a.c3, a.c4 + num)
Base.:(-)(num::Number,a::Even) = (-a.c1 + num, -a.c2, -a.c3, -a.c4 + num)

#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.c1, num*a.c2, num*a.c3, num*a.c4)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.c1, num*a.c2, num*a.c3, num*a.c4)

function Base.:(*)(a::Even,b::Even) 
    Even( a.c1*b.c1 + a.c2*b.c3,
            a.c1*b.c2 + a.c2*b.c4,
            a.c3*b.c1 + a.c4*b.c3,
            a.c4*b.c4 + a.c3*b.c2)
end
    
function Base.:(*)(a::Even,b::Odd)
    Odd( a.c1*b.c1 + a.c2*b.c3,
            a.c1*b.c2 + a.c2*b.c4,
            a.c3*b.c1 + a.c4*b.c3,
            a.c4*b.c4 + a.c3*b.c2)   
end


function Base.:(*)(a::Odd, b::Even)
    Odd(  a.c1*conj(b.c4) - a.c2*conj(b.c2),
            - a.c1*conj(b.c3) + a.c2*conj(b.c1),
            a.c3*conj(b.c4) - a.c4*conj(b.c2),
            a.c4*conj(b.c1) - a.c3*conj(b.c3))    
end


function Base.:(*)(a::Odd, b::Odd)
    Even(  a.c1*conj(b.c4) - a.c2*conj(b.c2),
            - a.c1*conj(b.c3) + a.c2*conj(b.c1),
            a.c3*conj(b.c4) - a.c4*conj(b.c2),
            a.c4*conj(b.c1) - a.c3*conj(b.c3))   
end


#Reverse
LinearAlgebra.adjoint(a::Even) = Even(a.c4, -a.c2, -a.c3, a.c1)
LinearAlgebra.adjoint(a::Odd) = Odd(conj(a.c1), conj(a.c3), conj(a.c2), conj(a.c4))


#Grade and projection
function project(a::Even,n::Integer)
    tra = (a.c1+a.c4)/2
    if (n==0)
        return Even((tra+conj(tra))/2, zero(a.c1), zero(a.c1), (tra+conj(tra))/2)
    elseif (n==2)
        return Even((a.c1-a.c4)/2, a.c2, a.c3, (a.c4-a.c1)/2 )
    elseif (n==4)
        return Even((tra-conj(tra))/2, zero(a.c1), zero(a.c1), (tra-conj(tra))/2)
    else
        return Even(zero(a.c1), zero(a.c1), zero(a.c1), zero(a.c1))
    end
end

function project(a::Odd,n::Integer)
    if (n==1)
        return Odd((a.c1+conj(a.c1))/2, (a.c2+conj(a.c3))/2, (a.c3+conj(a.c2))/2 , (a.c4+conj(a.c4))/2)
    elseif (n==3)
        return Odd((a.c1-conj(a.c1))/2, (a.c2-conj(a.c3))/2, (a.c3-conj(a.c2))/2 , (a.c4-conj(a.c4))/2)
    else
        return Odd(zero(a.c1), zero(a.c1), zero(a.c1), zero(a.c1))
    end
end


LinearAlgebra.tr(a::Even) = (real(a.c1+a.c4))/2
LinearAlgebra.dot(a::Even, b::Even) = real(a.c1*b.c1 + a.c2*b.c3 + a.c4*b.c4 + a.c3*b.c2)/2    
LinearAlgebra.dot(a::Odd, b::Odd) = real(a.c1*conj(b.c4) - a.c2*conj(b.c2) + a.c4*conj(b.c1) - a.c3*conj(b.c3))/2   


#Exponentiation
function expb(a::Even)
    a = project(a,2)
    aa = a*a
    fct = sqrt((aa.c1+aa.c4)/2)
    if iszero(fct)
        return 1+a
    else
        return cosh(fct)+sinh(fct)/fct*a
    end
end

function Base.exp(a::Even)
    R = expb(a)
    fct = (a.c1+a.c4)/2
    if iszero(fct)
        return R
    else
        return exp(fct)*R
    end
end

#Comparison. Uses default tolerances
Base.isapprox(a::Even, b::Even) = isapprox(a.c1,b.c1) && isapprox(a.c2,b.c2)  && isapprox(a.c3,b.c3) && isapprox(a.c4,b.c4) 
Base.isapprox(a::Odd, b::Odd) = isapprox(a.c1,b.c1) && isapprox(a.c2,b.c2)  && isapprox(a.c3,b.c3) && isapprox(a.c4,b.c4) 
