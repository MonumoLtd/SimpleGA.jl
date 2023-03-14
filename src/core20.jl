
#=
Core code for the implementation of GA(2,0)
Underlying representation is with complex numbers, so essentially a wrapper for Julia's internal ComplexF{T} formats.
=#

import ..project
import ..expb
import LinearAlgebra.tr
import LinearAlgebra.dot
import LinearAlgebra.adjoint


struct Even{T<:Real} <: Number
    c1::Complex{T}
end


struct Odd{T<:Real} <: Number
    c1::Complex{T}
end


#Addition / subtraction
Base.:(-)(a::Even) = Even(-a)
Base.:(-)(a::Odd) = Odd(-a)
Base.:(+)(a::Even,b::Even) = Even(a.c1 + b.c1)
Base.:(+)(a::Odd,b::Odd) = Odd(a.c1 + b.c1)
Base.:(-)(a::Even,b::Even) = Even(a.c1 - b.c1)
Base.:(-)(a::Odd,b::Odd) = Odd(a.c1 - b.c1)


#Scalar addition / subtraction. Other cases are in GAcommon
#Relies on Julia's promotion rules to do the sensible thing.
Base.:(+)(num::Number,a::Even) = Even(a.c1 + num)
Base.:(-)(num::Number,a::Even) = Even(-a.c1 + num)
   

#Multiplication
Base.:(*)(num::Number,a::Even) = Even(num*a.c1)
Base.:(*)(num::Number,a::Odd) = Odd(num*a.c1)
Base.:(*)(a::Even,b::Even) = Even( a.c1*b.c1) 
Base.:(*)(a::Even,b::Odd) = Odd(conj(a.c1)*b.c1)
Base.:(*)(a::Odd, b::Even) = Odd(a.c1*b.c1)   
Base.:(*)(a::Odd, b::Odd) = Even(conj(a.c1)*b.c1)    


#Division by a real
Base.:(/)(a::Even,num::Number) = Even((1/num)*a.c1)
Base.:(/)(a::Odd,num::Number) = Odd((1/num)*a.c1)

#Allow this here as complex numbers are a division algebra.
Base.:(/)(a::Even,b::Even) = Even(a.c1/b.c1)


#Reverse
adjoint(a::Even) = Even(conj(a.c1))
adjoint(a::Odd) = a

#Grade and projection
function project(a::Even,n::Integer)
    if (n==0)
        return Even(real(a.c1))
    elseif (n==2)
        return Even(imag(a.c1)*im )
    else
        return Even(zero(a.c1))
    end
end

function project(a::Odd,n::Integer)
    if (n==1)
        return a
    else
        return Odd(zero(a.c1))
    end
end

tr(a::Even) = real(a.c1)
dot(a::Even, b::Even) = real(a.c1*b.c1)    
dot(a::Odd, b::Odd) = real(conj(a.c1)*b.c1)


#Exponentiation
Base.exp(a::Even) = Even(exp(a.c1))

function expb(a::Even)
    Even(exp(im*a.c1.im))
end

#Comparison, using default tolerances.
Base.isapprox(a::Even, b::Even) = isapprox(a.c1,b.c1) 
Base.isapprox(a::Odd, b::Odd) = isapprox(a.c1,b.c1) 


