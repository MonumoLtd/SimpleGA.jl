#=
Code for GA(3,0). Even and odd elements are stored as quaternions.
=#

module GA30

include("GAcore30.jl")
include("GAcommon.jl")
import Base.show

const e1 = MVodd(0,1,0,0)
const e2 = MVodd(0,0,1,0)
const e3 = MVodd(0,0,0,1)
const I3 = MVodd(1,0,0,0)

bas30 = [e1,e2,e3]
export bas30

#Sets tolerance for not displaying results. Does not change the underlying multivector.
function approxzero(x::Float64)
    isapprox(x,0.0; atol = 1e-14)
end


function mvtype(a::MVeven)
    res=""
    tp = approxzero(a.w) ? "" : " + " * string(a.w)
    res *= tp
    tp = approxzero(a.x) ? "" : " + " * string(-a.x) * "e2e3"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(-a.y) * "e3e1"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(-a.z) * "e1e2"
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
    tp = approxzero(a.x) ? "" : " + " * string(a.x) * "e1"
    res *= tp
    tp = approxzero(a.y) ? "" : " + " * string(a.y) * "e2"
    res *= tp
    tp = approxzero(a.z) ? "" : " + " * string(a.z) * "e3"
    res *= tp
    tp = approxzero(a.w) ? "" : " + " * string(a.w) * "e123"
    res *= tp
    if (length(res) == 0)
        res = "0.0"
    else
        res = chop(res,head=3,tail=0)
    end
    return res
end

include("GAshow.jl")

end