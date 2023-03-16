#=
Code for GA(3,0,1). The Euclidean group algebra. 
=#

module PGA

include("PGAcore.jl")
include("GAcommon.jl")
import Base.show

#Sets tolerance for not displaying results. 
function approxzero(x::Float64)
    isapprox(x,0.0; atol = 1e-14)
end

#Basis
const e1 = MVodd(Quaternion(0, -1.0, 0, 0),qzero)
const e2 = MVodd(Quaternion(0, 0, -1.0, 0),qzero)
const e3 = MVodd(Quaternion(0, 0, 0, -1.0),qzero)
const e0 = MVodd(qzero, Quaternion(1.0,0,0,0))
const I3 = e1*e2*e3

basPGA = [e1,e2,e3, e0]
export basPGA, pdual

#Additional Functions. Sign conventions agree with De Keninck et al.
function pdual(a::MVeven)
    MVeven(-reverse(a.n),-reverse(a.q))
end

function pdual(a::MVodd)
    MVodd(-a.n,-a.q)
end

function mvtype(a::MVeven)
    res=""
    evenlist = [a.q.w, -a.q.z, - a.q.x, -a.q.y, -a.n.x, -a.n.y, -a.n.z, -a.n.w]
    evenstring = ["", "e1e2", "e2e3", "e3e1", "e0e1", "e0e2", "e0e3", "e0I3" ]
    for i = 1:8
        tp =  approxzero(evenlist[i]) ? "" : " + " * string(evenlist[i])*evenstring[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end


function mvtype(a::MVodd)
    res=""
    oddlist = [a.n.w,-a.q.x,-a.q.y,-a.q.z,a.n.z,a.n.y,a.n.x,-a.q.w]
    oddstring = ["e0", "e1", "e2", "e3", "e0e2e1", "e0e1e3", "e0e3e2", "I3"]
    for i = 1:8
        tp =  approxzero(oddlist[i]) ? "" : " + " * string(oddlist[i])*oddstring[i]
        res *= tp
    end
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

include("GAshow.jl")

end #Module