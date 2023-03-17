module GA24

#=
Code for G(2,4)
=#

using LinearAlgebra
using StaticArrays

include("core24.jl")
include("common.jl")

const g0 = Odd(id4)
const g1 = Odd(SMatrix{4,4,Complex{Int8},16}([0 1 0 0 ; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]))
const g2 = Odd(SMatrix{4,4,Complex{Int8},16}([0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]))
const g3 = Odd(SMatrix{4,4,Complex{Int8},16}([1 0 0 0; 0 -1 0 0 ; 0 0 -1 0; 0 0 0 1]))
const g4 = Odd(SMatrix{4,4,Complex{Int8},16}([0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]))
const g5 = Odd(SMatrix{4,4,Complex{Int8},16}([0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0 ]))

const I6 = g0*g1*g2*g3*g4*g5

bas24 = [1.0*g0, 1.0*g1, 1.0*g2, 1.0*g3, 1.0*g4, 1.0*g5]

function basis24(T) 
    scl = convert(T,1)
    return [scl*g0, scl*g1, scl*g2, scl*g3, scl*g4, scl*g5]
end

export bas24, basis24


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)

function mvtype(a::Even)
    scl = tr(a)
    res = approxzero(scl) ? "" : " + " * string(scl)
    tpevenbas = [g1*g0, g2*g0, g3*g0, g4*g0, -g5*g0, -g1*g2, -g1*g3, -g1*g4, g1*g5, -g2*g3, -g2*g4, g2*g5, -g3*g4, g3*g5, g4*g5,
        -g1*g0*I6, -g2*g0*I6, -g3*g0*I6, -g4*g0*I6, g5*g0*I6, g1*g2*I6, g1*g3*I6, g1*g4*I6, -g1*g5*I6, g2*g3*I6, g2*g4*I6, -g2*g5*I6, g3*g4*I6, -g3*g5*I6, -g4*g5*I6, -I6]
    tpevenstr = ["g1g0", "g2g0", "g3g0", "g4g0", "g5g0", "g1g2", "g1g3", "g1g4", "g1g5", "g2g3", "g2g4", "g2g5", "g3g4", "g3g5", "g4g5",
        "g1g0I6", "g2g0I6", "g3g0I6", "g4g0I6", "g5g0I6", "g1g2I6", "g1g3I6", "g1g4I6", "g1g5I6", "g2g3I6", "g2g4I6", "g2g5I6", "g3g4I6", "g3g5I6", "g4*g5*I6", "I6"]
    for i in 1:31
        scl = dot(a,tpevenbas[i])
        tp = approxzero(scl) ? "" : " + " * string(scl) * tpevenstr[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

function mvtype(a::Odd)
    tpoddbas = [g0, -g1, -g2, -g3, -g4, g5, 
        -g0*g1*g2, -g0*g1*g3, -g0*g1*g4, g0*g1*g5, -g0*g2*g3, -g0*g2*g4, g0*g2*g5, -g0*g3*g4, g0*g3*g5, g0*g4*g5,
        -g0*g1*g2*I6, -g0*g1*g3*I6, -g0*g1*g4*I6, g0*g1*g5*I6, -g0*g2*g3*I6, -g0*g2*g4*I6, g0*g2*g5*I6, -g0*g3*g4*I6, g0*g3*g5*I6, g0*g4*g5*I6,
        g0*I6, -g1*I6, -g2*I6, -g3*I6, -g4*I6, g5*I6 ]
    tpoddstr = ["g0", "g1", "g2", "g3", "g4", "g5", 
    "g0g1g2", "g0g1g3", "g0g1g4", "g0g1g5", "g0g2g3", "g0g2g4", "g0g2g5", "g0g3g4", "g0g3g5", "g0g4g5",
    "g0g1g2I6", "g0g1g3I6", "g0g1g4I6", "g0g1g5I6", "g0g2g3I6", "g0g2g4*I6", "g0g2g5I6", "g0g3g4I6", "g0g3g5I6", "g0g4g5I6",
    "g0I6", "g1I6", "g2I6", "g3I6", "g4I6", "g5I6"]
    res=""
    for i in 1:32
        scl = dot(a,tpoddbas[i])
        tp = approxzero(scl) ? "" : " + " * string(scl) * tpoddstr[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

include("show.jl")

end #Module


