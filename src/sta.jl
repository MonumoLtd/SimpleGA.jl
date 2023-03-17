
module STA
#=
Code for GA(1,3). Even and odd elements are stored as complex 'matrices'. 
Gamma and sigma used for pretty typing.
=#

using LinearAlgebra

include("coresta.jl")
include("common.jl")

#Basis
const s1 = Even{Int8}(0,1,1,0)
const s2 = Even{Int8}(0,-im,im,0)
const s3 = Even{Int8}(1,0,0,-1)
const g0 = Odd{Int8}(1,0,0,1)
const g1 = s1*g0
const g2 = s2*g0
const g3 = s3*g0
const I4 = g0*g1*g2*g3

basSTA = [g0,g1,g2,g3]

function basisSTA(T)
    s1 = Even{T}(0,1,1,0)
    s2 = Even{T}(0,-im,im,0)
    s3 = Even{T}(1,0,0,-1)
    g0 = Odd{T}(1,0,0,1)
    return [g0 , s1*g0, s2*g0, s3*g0 ]
end

export basSTA, basisSTA, bar

#Additional function for Pauli operations.Pre and post multiply by g0.
function bar(a::Even)
    Even(conj(a.c4), - conj(a.c3), - conj(a.c2), conj(a.c1))
end

function bar(a::Odd)
    Odd(conj(a.c4), - conj(a.c3), - conj(a.c2), conj(a.c1))
end


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)


function mvtype(a::Even)
    res=""
    scl = tr(a)
    tp = approxzero(scl) ? "" : " + " * string(scl)
    res *= tp
    scl = dot(a,s1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "σ1"
    res *= tp
    scl = dot(a,s2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "σ2"
    res *= tp
    scl = dot(a,s3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "σ3"
    res *= tp
    scl = dot(a,-I4*s1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iσ1"
    res *= tp
    scl = dot(a,-I4*s2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iσ2"
    res *= tp
    scl = dot(a,-I4*s3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iσ3"
    res *= tp
    scl = dot(a,-I4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end


function mvtype(a::Odd)
    res=""
    scl = dot(a,g0)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "γ0"
    res *= tp
    scl = dot(a,-g1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "γ1"
    res *= tp
    scl = dot(a,-g2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "γ2"
    res *= tp
    scl = dot(a,-g3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "γ3"
    res *= tp
    scl = dot(a,I4*g0)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iγ0"
    res *= tp
    scl = dot(a,-I4*g1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iγ1"
    res *= tp
    scl = dot(a,-I4*g2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iγ2"
    res *= tp
    scl = dot(a,-I4*g3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "Iγ3"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

include("show.jl")

end #Module