"""
    GA30

Module representing the `GA(3, 0)` geometric algebra.
"""
module GA30

using GADraft

struct MV{Parity,T} <: MultiVector{Parity,T}
    w::T
    x::T
    y::T
    z::T
end
function MV{P}(w, x, y, z) where {P<:Parity}
    w, x, y, z = promote(w, x, y, z)
    return MV{P,typeof(w)}(w, x, y, z)
end

# Default basis uses Float64.
const e1 = MV{Odd}(0.0, 1.0, 0.0, 0.0)
const e2 = MV{Odd}(0.0, 0.0, 1.0, 0.0)
const e3 = MV{Odd}(0.0, 0.0, 0.0, 1.0)
const I3 = MV{Odd}(1.0, 0.0, 0.0, 0.0)

end  # module GA30
