"""
    GA20

Module representing the `GA(2, 0)` geometric algebra.
"""
module GA20

# Export the default basis.
export e1, e2, I2

using GADraft
using GADraft: _complex_T, _scalar_T

struct MV{P,T} <: MultiVector{P,T}
    c1::Complex{T}
end
MV{P}(c1::Complex) where {P<:Parity} = MV{P,_complex_T(c1)}(c1)
MV{P}(r::Real) where {P<:Parity} = MV{P,typeof(r)}(complex(r))

# Construction
Base.zero(::MV{P,T}) where {P<:Parity,T} = MV{P}(zero(Complex{T}))

# Comparison
Base.isapprox(a::MV{P}, b::MV{P}; kwargs...) where {P} = isapprox(a.c1, b.c1; kwargs...)

# TODO Tom needs to understand what this means physically
# # Functions on even elements.
# Base.log(a::MV{Even}) = MV{Even}(log(a.c1))
# Base.real(a::MV{Even}) = MV{Even}(real(a.c1)) # not wrapped?
# Base.imag(a::MV{Even}) = MV{Even}(imag(a.c1)) # not wrapped?

function Base.show(io::IO, a::MV)
    res = if iseven(a)
        (iszero(real(a.c1)) ? "" : " + $(real(a.c1))") *
        (iszero(imag(a.c1)) ? "" : " + $(imag(a.c1)) I2")
    else  # isodd(a)
        (iszero(real(a.c1)) ? "" : " + $(real(a.c1)) e1") *
        (iszero(imag(a.c1)) ? "" : " + $(imag(a.c1)) e2")
    end
    res = isempty(res) ? "0.0" : res[4:end]
    return print(io, res)
end

# Unary negation.
Base.:(-)(a::MV{P}) where {P} = MV{P}(-a.c1)

# Addition & subtraction only defined between matching parities.
Base.:(+)(a::MV{P}, b::MV{P}) where {P} = MV{P}(a.c1 + b.c1)
Base.:(-)(a::MV{P}, b::MV{P}) where {P} = MV{P}(a.c1 - b.c1)

# Multiplication rules depend on parity.
Base.:(*)(a::MV{Even}, b::MV{Even}) = MV{Even}(a.c1 * b.c1)
Base.:(*)(a::MV{Even}, b::MV{Odd}) = MV{Odd}(conj(a.c1) * b.c1)
Base.:(*)(a::MV{Odd}, b::MV{Even}) = MV{Odd}(a.c1 * b.c1)
Base.:(*)(a::MV{Odd}, b::MV{Odd}) = MV{Even}(conj(a.c1) * b.c1)

# Addition and subtraction with scalars.
Base.:(+)(x::Real, a::MV{Even}) = MV{Even}(x + a.c1)
Base.:(+)(a::MV{Even}, x::Real) = MV{Even}(a.c1 + x)
Base.:(-)(x::Real, a::MV{Even}) = MV{Even}(x - a.c1)
Base.:(-)(a::MV{Even}, x::Real) = MV{Even}(a.c1 - x)

# Multiplication with scalars.
Base.:(*)(x::Real, a::MV{P}) where {P} = MV{P}(x * a.c1)
Base.:(*)(a::MV{P}, x::Real) where {P} = MV{P}(a.c1 * x)

"""
    reverse(a::MultiVector)

Reverse order of all vectors that might have appeared in a geometric product defining `a`.
"""
Base.reverse(a::MV{Even}) = MV{Even}(conj(a.c1))
Base.reverse(a::MV{Odd}) = a

function GADraft.project(a::MV{Even}, grade::Integer)
    return if grade == 0
        MV{Even}(scalar(a))
    elseif grade == 2
        MV{Even}(imag(a.c1) * im)
    else
        # Any other grade won't have a contribution.
        # TODO should we throw exceptions for non-sensical grades?
        #   (e.g. negative or larger than 2)
        zero(a)
    end

end

function GADraft.project(a::MV{Odd}, grade::Integer)
    return if grade == 1
        a
    else
        # Any other grade won't have a contribution.
        # TODO should we throw exceptions for non-sensical grades?
        zero(a)
    end
end

GADraft.scalar(a::MV{Even}) = real(a.c1)
GADraft.scalar(a::MV{Odd}) = zero(_scalar_T(a))
GADraft.scalar(a::MV{Even}, b::MV{Even}) = real(a.c1 * b.c1)
GADraft.scalar(a::MV{Odd}, b::MV{Odd}) = real(conj(a.c1) * b.c1)
# Catches Even * Odd and Odd * Even.
GADraft.scalar(a::MV, b::MV) = zero(promote_type(_scalar_T(a), _scalar_T(b)))

# Default basis uses Float64.
const e1 = MV{Odd,Float64}(1.0)
const e2 = MV{Odd,Float64}(im)
const I2 = MV{Even,Float64}(im)

end  # module GA20
