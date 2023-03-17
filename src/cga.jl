module CGA

#=
Code for GA(4,1). 
=#


using LinearAlgebra

include("corecga.jl")
include("common.jl")
import Base.show


#Basis
const e1 = Odd{Int8}(Quaternion(0,1,0,0),Quaternion(0.0),Quaternion(0.0),Quaternion(0,-1,0,0))
const e2 = Odd{Int8}(Quaternion(0,0,1,0),Quaternion(0.0),Quaternion(0.0),Quaternion(0,0,-1,0))
const e3 = Odd{Int8}(Quaternion(0,0,0,1),Quaternion(0.0),Quaternion(0.0),Quaternion(0,0,0,-1))
const e4 = Odd{Int8}(Quaternion(0.0),Quaternion(1,0,0,0),Quaternion(1,0,0,0),Quaternion(0.0))
const f4 = Odd{Int8}(Quaternion(0.0),Quaternion(-1,0,0,0),Quaternion(1,0,0,0),Quaternion(0.0))
const I5 = e1*e2*e3*e4*f4

function basisCGA(T)
    e1 = Odd{T}(Quaternion{T}(0,1,0,0),Quaternion{T}(0),Quaternion{T}(0),Quaternion{T}(0,-1,0,0))
    e2 = Odd{T}(Quaternion{T}(0,0,1,0),Quaternion{T}(0),Quaternion{T}(0),Quaternion{T}(0,0,-1,0))
    e3 = Odd{T}(Quaternion{T}(0,0,0,1),Quaternion{T}(0),Quaternion{T}(0),Quaternion{T}(0,0,0,-1))
    e4 = Odd{T}(Quaternion{T}(0),Quaternion{T}(1,0,0,0),Quaternion{T}(1,0,0,0),Quaternion{T}(0))
    f4 = Odd{T}(Quaternion{T}(0),Quaternion{T}(-1,0,0,0),Quaternion{T}(1,0,0,0),Quaternion{T}(0))
    return [e1,e2,e3,e4,f4]
end

basCGA = [e1,e2,e3,e4,f4]
export basCGA, basisCGA


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)


function mvtype(a::Even)
    res=""
    scl = tr(a)
    tp = approxzero(scl) ? "" : " + " * string(scl)
    res *= tp
    scl = dot(a,-e1*e2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e2"
    res *= tp
    scl = dot(a,-e2*e3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2e3"
    res *= tp
    scl = dot(a,-e3*e1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3e1"
    res *= tp
    scl = dot(a,-e1*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e4"
    res *= tp
    scl = dot(a,-e2*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2e4"
    res *= tp
    scl = dot(a,-e3*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3e4"
    res *= tp
    scl = dot(a,e1*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1f4"
    res *= tp
    scl = dot(a,e2*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2f4"
    res *= tp
    scl = dot(a,e3*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3f4"
    res *= tp
    scl = dot(a,e4*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e4f4"
    res *= tp
    scl = dot(a,-e1*I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5e1"
    res *= tp
    scl = dot(a,-e2*I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5e2"
    res *= tp
    scl = dot(a,-e3*I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5e3"
    res *= tp
    scl = dot(a,-e4*I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5e4"
    res *= tp
    scl = dot(a,f4*I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5f4"
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
    scl = dot(a,e1)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1"
    res *= tp
    scl = dot(a,e2)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2"
    res *= tp
    scl = dot(a,e3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3"
    res *= tp
    scl = dot(a,e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e4"
    res *= tp
    scl = dot(a,-f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "f4"
    res *= tp
    scl = dot(a,-e1*e2*e3)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e2e3"
    res *= tp
    scl = dot(a,-e1*e2*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e2e4"
    res *= tp
    scl = dot(a,-e1*e3*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e3e4"
    res *= tp
    scl = dot(a,-e2*e3*e4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2e3e4"
    res *= tp
    scl = dot(a,e1*e2*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e2f4"
    res *= tp
    scl = dot(a,e1*e3*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e3f4"
    res *= tp
    scl = dot(a,e2*e3*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2e3f4"
    res *= tp
    scl = dot(a,e1*e4*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1e4f4"
    res *= tp
    scl = dot(a,e2*e4*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2e4f4"
    res *= tp
    scl = dot(a,e3*e4*f4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3e4f4"
    res *= tp
    scl = dot(a,-I5)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "I5"
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