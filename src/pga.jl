
module PGA

#=
Code for GA(3,0,1). The Euclidean group algebra. 
=#


using LinearAlgebra

include("corepga.jl")
include("common.jl")


#Basis
const e1 = Odd{Float64}(Quaternion(0, -1.0, 0, 0), Quaternion(0.0))
const e2 = Odd{Float64}(Quaternion(0, 0, -1.0, 0),Quaternion(0.0))
const e3 = Odd{Float64}(Quaternion(0, 0, 0, -1.0),Quaternion(0.0))
const e0 = Odd{Float64}(Quaternion(0.0), Quaternion(1.0,0,0,0))
const I3 = e1*e2*e3

basPGA = [e1,e2,e3, e0]

function basisPGA(T)
    e1 = Odd{T}(Quaternion{T}(0, -1, 0, 0), Quaternion{T}(0))
    e2 = Odd{T}(Quaternion{T}(0, 0, -1, 0),Quaternion{T}(0))
    e3 = Odd{T}(Quaternion{T}(0, 0, 0, -1),Quaternion{T}(0))
    e0 = Odd{T}(Quaternion{T}(0), Quaternion{T}(1,0,0,0))
    return [e1,e2,e3,e0]
end

export basPGA, basisPGA, pdual


#Additional Functions. Sign conventions agree with De Keninck et al.
function pdual(a::Even)
    Even(-reverse(a.n),-reverse(a.q))
end

function pdual(a::Odd)
    Odd(-a.n,-a.q)
end


#Sets tolerance for not displaying results. Adding 1 to comparison seems to work well.
approxzero(x::Real) = isapprox(1+x,1.0)

function mvtype(a::Even)
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


function mvtype(a::Odd)
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

include("show.jl")

end #Module