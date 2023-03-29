"""
    GA30

Module representing the GA(3, 0) geometric algebra.
Underlying rep is quaternions, but not used explicitly here to keep this code self-contained.
"""

module GA30

using GeometricAlgebra
using LinearAlgebra
using StaticArrays

include("core30.jl")
include("common.jl")

#Basis
const e1 = Odd{Int8}(0, 1, 0, 0)
const e2 = Odd{Int8}(0, 0, 1, 0)
const e3 = Odd{Int8}(0, 0, 0, 1)
const I3 = Odd{Int8}(1, 0, 0, 0)

const basis = SA[e1, e2, e3]

#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1 + x, 1.0)

function mv_to_text(a::Even)
    res = ""
    tp = approxzero(a.w) ? "" : " + " * string(a.w)
    res *= tp
    tp = approxzero(a.x) ? "" : " + " * string(-a.x) * "e2e3"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(-a.y) * "e3e1"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(-a.z) * "e1e2"
    res *= tp
    if (length(res) == 0)
        res = string(zero(a.w))
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

function mv_to_text(a::Odd)
    res = ""
    tp = approxzero(a.x) ? "" : " + " * string(a.x) * "e1"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(a.y) * "e2"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(a.z) * "e3"
    res *= tp
    tp = approxzero(a.w) ? "" : " + " * string(a.w) * "e123"
    res *= tp
    if (length(res) == 0)
        res = string(zero(a.w))
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

include("show.jl")

end
