"""
    GA20

Module representing the `GA(2, 0)` geometric algebra.
Even and odd elements are stored as complex numbers.

"""
module GA20
using LinearAlgebra

include("core20.jl")
include("common.jl")

#Basis
const e1 = Odd{Int8}(1)
const e2 = Odd{Int8}(im)
const I2 = Even{Int8}(im)

const basis = [e1, e2]

#Extra functions on even elements
Base.log(a::Even) = Even(log(a.c1))
Base.real(a::Even) = real(a.c1)
Base.imag(a::Even) = imag(a.c1)

#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1 + x, 1.0)

function mv_to_text(a::Even)
    res = ""
    tp = approxzero(real(a.c1)) ? "" : " + " * string(real(a.c1))
    res *= tp
    tp = approxzero(imag(a.c1)) ? "" : " + " * string(imag(a.c1)) * "I2"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

function mv_to_text(a::Odd)
    res = ""
    tp = approxzero(real(a.c1)) ? "" : " + " * string(real(a.c1)) * "e1"
    res *= tp
    tp = approxzero(imag(a.c1)) ? "" : " + " * string(imag(a.c1)) * "e2"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

include("show.jl")

end  #Module