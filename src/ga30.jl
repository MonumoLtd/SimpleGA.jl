"""
    GA30

Module representing the GA(3, 0) geometric algebra.
Underlying rep is quaternions, but not used explicitly here to keep this code self-contained.
"""
module GA30

using SimpleGA
using LinearAlgebra
using StaticArrays

include("core30.jl")
include("common.jl")

# Basis
const e1 = Odd{Int8}(0, 1, 0, 0)
const e2 = Odd{Int8}(0, 0, 1, 0)
const e3 = Odd{Int8}(0, 0, 0, 1)
const I3 = Odd{Int8}(1, 0, 0, 0)

const basis = SA[e1, e2, e3]

function mv_to_text(a::Even)
    res = ""
    tp = iszero(a.w) ? "" : " + " * string(a.w)
    res *= tp
    tp = iszero(a.x) ? "" : " + " * string(-a.x) * "e2e3"
    res *= tp
    tp = iszero(a.y) ? "" : " + " * string(-a.y) * "e3e1"
    res *= tp
    tp = iszero(a.z) ? "" : " + " * string(-a.z) * "e1e2"
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
    tp = iszero(a.x) ? "" : " + " * string(a.x) * "e1"
    res *= tp
    tp = iszero(a.y) ? "" : " + " * string(a.y) * "e2"
    res *= tp
    tp = iszero(a.z) ? "" : " + " * string(a.z) * "e3"
    res *= tp
    tp = iszero(a.w) ? "" : " + " * string(a.w) * "e123"
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
