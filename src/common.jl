#=
All common code shared between the specific implementations is placed here.
Mostly this is code to control multiple dispatch for the various ways of combining with scalars.
=#

#Addition. Can only add a scalar to an even multivector.
Base.:(+)(a::Even,num::Number) = num + a
Base.:(-)(a::Even,num::Number) = (-num) + a

#Multiplication
Base.:(*)(a::Even,num::Number) = num * a
Base.:(*)(a::Odd,num::Number) = num * a

#Projection
tr(a::Odd) = 0
dot(a::Even,b::Odd) = 0
dot(a::Odd,b::Even) = 0