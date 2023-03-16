
module GA40

#=
Code for GA(4,0). Even and odd elements are stored as quaternion pairs.
=#

using LinearAlgebra

include("core40.jl")
include("common.jl")


#Basis
const e1 = Odd{Float64}(Quaternion(0,1,0,0), Quaternion(0,-1,0,0))
const e2 = Odd{Float64}(Quaternion(0,0,1,0), Quaternion(0,0,-1,0))
const e3 = Odd{Float64}(Quaternion(0,0,0,1), Quaternion(0,0,0,-1))
const e4 = Odd{Float64}(Quaternion(1,0,0,0), Quaternion(1,0,0,0))
const E4 = e1*e2*e3*e4

bas40 = [e1,e2,e3,e4]

function basis40(T)
    e1 = Odd{T}(Quaternion(0,1,0,0), Quaternion(0,-1,0,0))
    e2 = Odd{T}(Quaternion(0,0,1,0), Quaternion(0,0,-1,0))
    e3 = Odd{T}(Quaternion(0,0,0,1), Quaternion(0,0,0,-1))
    e4 = Odd{T}(Quaternion(1,0,0,0), Quaternion(1,0,0,0))
    return [e1,e2,e3,e4]
end
    
export bas40, basis40


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
    scl = dot(a,E4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "E4"
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
    scl = dot(a,-e1*E4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e1E4"
    res *= tp
    scl = dot(a,-e2*E4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e2E4"
    res *= tp
    scl = dot(a,-e3*E4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e3E4"
    res *= tp
    scl = dot(a,-e4*E4)
    tp = approxzero(scl) ? "" : " + " * string(scl) * "e4E4"
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