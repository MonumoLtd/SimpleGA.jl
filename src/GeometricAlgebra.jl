module GeometricAlgebra

using LinearAlgebra

# Prototype functions.
#Definitions are provided in each of the submodules
function project() end
function bivector_exp() end

export project, bivector_exp, inject

#Useful function for doting vectors into GA bases.
#Dot works in most cases, but not for GA44.
#TODO: Understand why dot() does not work for GA44.
function inject(xs, ys)
    return reduce(+, map((x, y) -> x * y, xs, ys))
end

#Each algebra is in a separate sub-module.

include("quaternions.jl")
using .Quaternions

include("ga20.jl")
using .GA20
export GA20

include("ga30.jl")
using .GA30
export GA30

include("sta.jl")
using .STA
export bar, STA

include("ga44.jl")
using .GA44
export construct44, GA44

include("ga40.jl")
using .GA40
export GA40

include("pga.jl")
using .PGA
export pdual, PGA

include("cga.jl")
using .CGA
export CGA

include("ga33.jl")
using .GA33
export GA33

include("ga24.jl")
using .GA24
export GA24

include("ga3232.jl")
using .GA3232
export construct64, GA3232

# The main basis function.
include("basis.jl")
export basis, testbas

# Maps into GA(4,4). Primarily for testing.
include("embed.jl")
export embed

end #module
