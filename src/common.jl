#=
All common code shared between the specific implementations is placed here.
Mostly this is code to control multiple dispatch for the various ways of combining with scalars.
=#

#Addition. Can only add a scalar to an even multivector.
Base.:(+)(a::Even, num::Number) = num + a
Base.:(-)(a::Even, num::Number) = (-num) + a

#Multiplication
Base.:(*)(a::Even, num::Number) = num * a
Base.:(*)(a::Odd, num::Number) = num * a

#Division by a real
Base.:(/)(a::Even, num::Number) = (1 / num) * a
Base.:(/)(a::Odd, num::Number) = (1 / num) * a

#Projection
LinearAlgebra.tr(a::Odd) = 0
LinearAlgebra.dot(a::Even, b::Odd) = 0
LinearAlgebra.dot(a::Odd, b::Even) = 0
