
using LinearAlgebra

module GeometricAlgebra


# Core functions.
export project, expb, inject, basis

function project() end
function expb() end

function inject(xs,ys) 
    reduce(+,map((x,y)->x*y,xs,ys))
end


#Each algebra in a separate sub-module.

include("quaternion.jl")
using .Quaternions

include("ga20.jl")
using .GA20

include("GA30.jl")
using .GA30

include("sta.jl")
using .STA
export bar

#=


include("GA40.jl")
using .GA40

include("PGA.jl")
using .PGA
export pdual

include("CGA.jl")
using .CGA


include("GA33.jl")
using .GA33


include("GA24.jl")
using .GA24

include("GA44.jl")
using .GA44
export construct44

include("GA64.jl")
using .GA64
export construct64

=#

function basis(alg)
    if alg == "GA20" 
        return bas20
    elseif alg =="GA30"
        return bas30
    elseif alg =="GA40"
        return bas40
    elseif alg =="PGA"
        return basPGA
    elseif alg =="CGA"
        return basCGA
    elseif alg =="STA"
        return basSTA
    elseif alg =="GA33"
        return bas33
    elseif alg == "GA24"
        return bas24
    elseif alg =="GA44"
        return bas44
    elseif alg == "GA64"
        return bas64
    end
end

function basis(T,alg)
    if alg == "GA20" 
        return basis20(T)
    elseif alg =="GA30"
        return basis30(T)
    elseif alg =="STA"
        return basisSTA(T)
    end
end

#=
include("embed.jl")
export embed
=#

end #module