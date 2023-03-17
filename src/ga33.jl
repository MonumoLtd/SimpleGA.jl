module GA33

#=
Code for GA(3,3). 
=#

using LinearAlgebra
using StaticArrays

include("core33.jl")
include("common.jl")


#Basis
const e1 = Odd(
    SMatrix{4,4,Int8}([0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]),
    SMatrix{4,4,Int8}([0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0])
)
const e2 = Odd(
    SMatrix{4,4,Int8}([0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0]),
    SMatrix{4,4,Int8}([0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0])
)
const e3 = Odd(
    id4, id4
)
const f1 = Odd(
    SMatrix{4,4,Int8}([0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0]),
    SMatrix{4,4,Int8}([0 0 -1 0; 0 0 0 1; -1 0 0 0; 0 1 0 0])
)
const f2 = Odd(
    SMatrix{4,4,Int8}([0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]),
    SMatrix{4,4,Int8}([0 0 0 -1; 0 0 -1 0; 0 -1 0 0; -1 0 0 0])
)
const f3 = Odd(
    SMatrix{4,4,Int8}([1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]),
    SMatrix{4,4,Int8}([-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1])
)

E6 = f1*e1*f2*e2*f3*e3

bas33 = [e1,e2,e3,f1,f2,f3]

function basis33(T)
    e1 = Odd(SMatrix{4,4,T}([0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]), SMatrix{4,4,T}([0 0 -1 0; 0 0 0 -1; 1 0 0 0; 0 1 0 0]))
    e2 = Odd(SMatrix{4,4,T}([0 0 0 1; 0 0 -1 0; 0 1 0 0; -1 0 0 0]), SMatrix{4,4,T}([0 0 0 -1; 0 0 1 0; 0 -1 0 0; 1 0 0 0]))
    e3 = Odd(id4(T), id4(T))
    f1 = Odd(SMatrix{4,4,T}([0 0 1 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0]),SMatrix{4,4,T}([0 0 -1 0; 0 0 0 1; -1 0 0 0; 0 1 0 0]))
    f2 = Odd(SMatrix{4,4,T}([0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]),SMatrix{4,4,T}([0 0 0 -1; 0 0 -1 0; 0 -1 0 0; -1 0 0 0]))
    f3 = Odd(SMatrix{4,4,T}([1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]),SMatrix{4,4,T}([-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]))
    return [e1,e2,e3,f1,f2,f3]
end

export bas33, basis33


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)


#Arrays to simplify MVtype code

function mvtype(a::Even)
    scl = tr(a)
    res = approxzero(scl) ? "" : " + " * string(scl)
    tpevenbas = [-e1*e2, -e1*e3, -e2*e3, -f1*f2, -f1*f3, -f2*f3, e1*f1, e1*f2, e1*f3, e2*f1, e2*f2, e2*f3, e3*f1, e3*f2, e3*f3,
    -e1*e2*E6, -e1*e3*E6, -e2*e3*E6, -f1*f2*E6, -f1*f3*E6, -f2*f3*E6, e1*f1*E6, e1*f2*E6, e1*f3*E6,e2*f1*E6, e2*f2*E6, e2*f3*E6, e3*f1*E6, e3*f2*E6, e3*f3*E6,E6]
    tpevenstr = ["e1e2", "e1e3", "e2e3", "f1f2", "f1f3", "f2f3", "e1f1", "e1f2", "e1f3", "e2f1", "e2f2", "e2f3", "e3f1", "e3f2", "e3f3",
    "e1e2E6", "e1e3E6", "e2e3E6", "f1f2E6", "f1f3E6", "f2f3E6", "e1f1E6", "e1f2E6", "e1f3E6", "e2f1E6", "e2f2E6", "e2f3E6", "e3f1E6", "e3f2E6", "e3f3E6","E6"]
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
    tpoddbas = [e1, e2, e3, -f1, -f2, -f3, -e1*e2*e3, 
    e1*e2*f1, e1*e3*f1, e2*e3*f1, e1*e2*f2, e1*e3*f2, e2*e3*f2, e1*e2*f3, e1*e3*f3, e2*e3*f3,
    -e1*f1*f2, - e1*f1*f3, -e1*f2*f3, -e2*f1*f2, - e2*f1*f3, -e2*f2*f3, -e3*f1*f2, - e3*f1*f3, -e3*f2*f3,
    f1*f2*f3, -e1*E6, -e2*E6, -e3*E6, f1*E6, f2*E6, f3*E6]
    tpoddstr = ["e1", "e2", "e3", "f1", "f2", "f3", "e1e2e3", 
    "e1e2f1", "e1e3f1", "e2e3f1", "e1e2f2", "e1e3f2", "e2e3f2", "e1e2f3", "e1e3f3", "e2e3f3",
    "e1f1f2", " e1f1f3", "e1f2f3", "e2f1f2", " e2f1f3", "e2f2f3", "e3f1f2", " e3f1f3", "e3f2f3",
    "f1f2f3", "e1E6", "e2E6", "e3E6", "f1E6", "f2E6", "f3E6"]
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