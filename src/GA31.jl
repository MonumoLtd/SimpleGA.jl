"""
    GA31

Code for GA(3,1). Even and odd elements are stored as complex structs.
2x2 matrix products are unwrapped.
"""
 
module GA31

using GeometricAlgebra
using LinearAlgebra
using StaticArrays

include("core31.jl")
include("common.jl")

#Basis
const s1 = Even{Int8}(0, 1, 1, 0)
const s2 = Even{Int8}(0, -im, im, 0)
const s3 = Even{Int8}(1, 0, 0, -1)
const f3 = Odd{Int8}(1, 0, 0, 1)
const e1 = s1 * f3
const e2 = s2 * f3
const e3 = s3 * f3
const I4 = e1 * e2 * e3 * f3

const basis = SA[e1, e2, e3, f3]

#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1 + x, 1.0)

function mv_to_text(a::Even)
    scl = tr(a)
    res = approxzero(scl) ? "" : " + " * string(scl)
     #! format:off
    tpevenbas = [-e1*e2, -e1*e3, -e2*e3, e1*f3 , e2*f3, e3*f3,
                -I4]
    tpevenstr = ["e1e2", "e1e3", "e2e3", "e1f3" , "e2f3", "e3f3",
                    "I4"]
    #! format:on
    for i in 1:7
        scl = dot(a, tpevenbas[i])
        tp = approxzero(scl) ? "" : " + " * string(scl) * tpevenstr[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

function mv_to_text(a::Odd)
    #! format:off
    tpoddbas = [e1, e2, e3, -f3,
                I4*e1, I4*e2, I4*e3, -I4*f3]
    tpoddstr = ["e1", "e2", "e3", "f3", 
                "I4e1", "I4e2", "I4e3", "I4f3"]
    #! format:on
    res = ""
    for i in 1:8
        scl = dot(a, tpoddbas[i])
        tp = approxzero(scl) ? "" : " + " * string(scl) * tpoddstr[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res; head=3, tail=0)
    end
    return res
end

include("show.jl")

end #Module