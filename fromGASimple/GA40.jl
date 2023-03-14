#=
Code for GA(4,0). Even and odd elements are stored as quaternion pairs.
=#

module GA40

include("GAcore40.jl")
include("GAcommon.jl")
import Base.show

#Sets tolerance for not displaying results. 
function approxzero(x::Float64)
    isapprox(x,0.0; atol = 1e-14)
end

#Basis
const e1 = MVodd(Quaternion(0,1,0,0), Quaternion(0,-1,0,0))
const e2 = MVodd(Quaternion(0,0,1,0), Quaternion(0,0,-1,0))
const e3 = MVodd(Quaternion(0,0,0,1), Quaternion(0,0,0,-1))
const e4 = MVodd(Quaternion(1,0,0,0), Quaternion(1,0,0,0))
const E4 = e1*e2*e3*e4

bas40 = [e1,e2,e3,e4]
export bas40

function mvtype(a::MVeven)
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


function mvtype(a::MVodd)
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

include("GAshow.jl")

end #Module