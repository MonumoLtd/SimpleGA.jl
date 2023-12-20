# All common code shared between the specific implementations is placed here.
# By convention it will be included into each GA-specific submodule

# Addition. Can only add a scalar to an even multivector.
Base.:(+)(a::Even, num::Real) = num + a
Base.:(-)(a::Even, num::Real) = (-num) + a

# Multiplication
Base.:(*)(a::Even, num::Real) = num * a
Base.:(*)(a::Odd, num::Real) = num * a

# Division by a real
Base.:(/)(a::Even, num::Real) = (1 / num) * a
Base.:(/)(a::Odd, num::Real) = (1 / num) * a

# Projection
LinearAlgebra.tr(a::Odd) = 0
LinearAlgebra.dot(a::Even, b::Odd) = 0
LinearAlgebra.dot(a::Odd, b::Even) = 0

# Approximate equality.
Base.isapprox(::Odd, ::Even; kwargs...) = false
Base.isapprox(::Even, ::Odd; kwargs...) = false
Base.isapprox(::AbstractArray{<:Odd}, ::AbstractArray{<:Even}; kwargs...) = false
Base.isapprox(::AbstractArray{<:Even}, ::AbstractArray{<:Odd}; kwargs...) = false

function Base.isapprox(a::AbstractArray{<:Even}, b::AbstractArray{<:Even}; kwargs...)
    size(a) != size(b) && return false
    return all(isapprox.(a, b; kwargs...))
end
function Base.isapprox(a::AbstractArray{<:Odd}, b::AbstractArray{<:Odd}; kwargs...)
    size(a) != size(b) && return false
    return all(isapprox.(a, b; kwargs...))
end
