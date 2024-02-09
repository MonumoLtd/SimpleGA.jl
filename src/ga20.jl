"""
    GA20

Module representing the `GA(2, 0)` geometric algebra.
Even and odd elements are stored as complex numbers.
"""
module GA20

using SimpleGA
using LinearAlgebra
using StaticArrays

include("core20.jl")
include("common.jl")

#Basis
const e1 = Odd{Int8}(1)
const e2 = Odd{Int8}(im)
const I2 = Even{Int8}(im)

const basis = SA[e1, e2]

#Extra functions on even elements
Base.log(a::Even) = Even(log(a.c1))
Base.real(a::Even) = real(a.c1)
Base.imag(a::Even) = imag(a.c1)

function mv_to_text(a::Even)
    res = ""
    tp = iszero(real(a.c1)) ? "" : " + " * string(real(a.c1))
    res *= tp
    tp = iszero(imag(a.c1)) ? "" : " + " * string(imag(a.c1)) * "I2"
    res *= tp
    if (length(res) == 0)
        res = string(zero(real(a.c1)))
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

function mv_to_text(a::Odd)
    res = ""
    tp = iszero(real(a.c1)) ? "" : " + " * string(real(a.c1)) * "e1"
    res *= tp
    tp = iszero(imag(a.c1)) ? "" : " + " * string(imag(a.c1)) * "e2"
    res *= tp
    if (length(res) == 0)
        res = string(zero(real(a.c1)))
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

include("show.jl")

end  #Module
