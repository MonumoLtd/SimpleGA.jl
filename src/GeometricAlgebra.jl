using LinearAlgebra

module GeometricAlgebra

# Utility functions.

function project() end
function expb() end

function inject(xs, ys)
    return reduce(+, map((x, y) -> x * y, xs, ys))
end

export project, expb, inject

#Each algebra is in a separate sub-module.

include("quaternion.jl")
using .Quaternions

include("ga20.jl")
using .GA20
export bas20

include("GA30.jl")
using .GA30
export bas30

include("sta.jl")
using .STA
export bar
export basSTA

include("ga44.jl")
using .GA44
export construct44, bas44

include("ga40.jl")
using .GA40
export bas40

include("PGA.jl")
using .PGA
export pdual, basPGA

include("CGA.jl")
using .CGA
export basCGA

include("GA33.jl")
using .GA33
export bas33

include("GA24.jl")
using .GA24
export bas24

include("GA64.jl")
using .GA64
export construct64, bas64

#The main basis function.
include("basis.jl")
export basis, testbas

#Maps into GA(4,4). Primarily for testing.
include("embed.jl")
export embed

end #module