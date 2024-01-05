#=
Implementation of GA(32,32) in Julia using Uint64 bitwise operations.
=#

#The Multivector type assumes that the blade list is unique and in order. But we want to avoid checking this at runtime.
#Only use this constructor if you are certain the blade list is correct. If not, use construct64()
struct Multivector{T<:Real} <: Number
    bas::Vector{UInt64}
    val::Vector{T}
end

# NOTE: we do not define equality in terms of the individual components, since we have many ways to represent "zero" that
#   we wish to be equivalent.
Base.:(==)(mv1::Multivector, mv2::Multivector) = iszero((mv1 - mv2).val)
Base.hash(a::Multivector, h::UInt) = hash(a.bas, hash(a.val, hash(:Multivector, h)))

function construct64(bs, vs)
    if length(bs) != length(unique(bs))
        error("List of blades must be unique")
    end

    p = sortperm(bs)
    return Multivector(bs[p], vs[p])
end

basscl = 0x0000000000000000

function Base.convert(::Type{Multivector{T}}, a::Multivector) where {T<:Real}
    return Multivector{T}(a.bas, convert.(T, a.val))
end
function Base.promote_rule(
    ::Type{Multivector{S}}, ::Type{Multivector{T}}
) where {S<:Real,T<:Real}
    return Multivector{promote_type(S, T)}
end

Base.zero(::Type{Multivector{T}}) where {T} = Multivector([basscl], [zero(T)])
Base.one(::Type{Multivector{T}}) where {T} = Multivector([basscl], [one(T)])
Base.zero(a::Multivector) = zero(typeof(a))
Base.one(a::Multivector) = one(typeof(a))

#This avoids rounding error bloating multivectors. Set of FP64
const mvtol = 1e-14

function mvtidy(mv::Multivector)
    ln = length(filter(x -> !isapprox(x, 0.0; atol=mvtol), mv.val))
    iszero(ln) && return zero(mv)

    rsbas = zeros(UInt64, ln)
    rsval = zeros(typeof(mv.val[1]), ln)
    j = 1
    for i in 1:length(mv.bas)
        if !isapprox(mv.val[i], 0.0)
            rsbas[j] = mv.bas[i]
            rsval[j] = mv.val[i]
            j += 1
        end
    end
    return Multivector(rsbas, rsval)
end

#Addition / subtraction
Base.:(-)(mv::Multivector) = Multivector(mv.bas, -mv.val)

function Base.:(+)(mv1::Multivector, mv2::Multivector)
    rsbas = sort(union(mv1.bas, mv2.bas))
    l1 = length(mv1.bas)
    l2 = length(mv2.bas)
    ln = length(rsbas)
    rsval = zeros(typeof(mv1.val[1] + mv2.val[1]), ln)
    i = 1
    j = 1
    for k in 1:ln
        if i <= l1 && rsbas[k] == mv1.bas[i]
            rsval[k] += mv1.val[i]
            i += 1
        end
        if j <= l2 && rsbas[k] == mv2.bas[j]
            rsval[k] += mv2.val[j]
            j += 1
        end
    end
    return Multivector(rsbas, rsval)
end

Base.:(+)(nm::Real, mv::Multivector) = mv + nm * one(mv)
Base.:(+)(mv::Multivector, nm::Real) = nm + mv
Base.:(-)(nm::Real, mv::Multivector) = nm + (-mv)
Base.:(-)(mv::Multivector, nm::Real) = (-nm) + mv
Base.:(-)(mv1::Multivector, mv2::Multivector) = mv1 + (-(mv2))

#Multiplication

Base.:(*)(num::Real, mv::Multivector) = Multivector(mv.bas, num * mv.val)
Base.:(*)(mv::Multivector, num::Real) = num * mv

function gaprodsign(bld1, bld2)
    tp1 = xor(bld2, bld2 << 1)
    cntones = count_ones((bld1 & 0xaaaaaaaaaaaaaaaa) & tp1)
    return convert(Int8, 1 - 2 * (cntones % 2))
end

function Base.:(*)(mv1::Multivector, mv2::Multivector)
    res = 0
    l1 = length(mv1.bas)
    for i in 1:l1
        bs = map(bd -> xor(mv1.bas[i], bd), mv2.bas)
        vs = map((bd, vl) -> mv1.val[i] * vl * gaprodsign(mv1.bas[i], bd), mv2.bas, mv2.val)
        p = sortperm(bs)
        res += Multivector(sort(bs), vs[p])
    end
    return mvtidy(res)
end

#Division by a real
Base.:(/)(mv::Multivector, num::Real) = (1 / num) * mv

#Reverse

#Converts conformal split to traditional rep for grade operation
function bldconvert(xin::UInt64)
    res = xin
    msk = 0xc000000000000000
    chk = ~msk
    for i in 1:31
        if count_ones(res & msk) == 1
            res = xor(res, chk)
        end
        msk = msk >> 2
        chk = chk >> 2
    end
    return res
end

grd(xin::UInt64) = count_ones(bldconvert(xin))

function LinearAlgebra.adjoint(mv::Multivector)
    rsval = similar(mv.val)
    for i in 1:length(mv.bas)
        if isodd(div(grd(mv.bas[i]), 2))
            rsval[i] = -mv.val[i]
        else
            rsval[i] = mv.val[i]
        end
    end
    return Multivector(mv.bas, rsval)
end

#Grade and projection
function SimpleGA.project(mv::Multivector, n::Int64)
    rsbas = filter(x -> grd(x) == n, mv.bas)
    ln = length(rsbas)
    iszero(ln) && return zero(mv)

    rsval = zeros(typeof(mv.val[1]), ln)
    for i in 1:ln
        j = findfirst(isequal(rsbas[i]), mv.bas)
        rsval[i] = mv.val[j]
    end
    return Multivector(rsbas, rsval)
end

LinearAlgebra.tr(mv::Multivector) = mv.bas[1] == basscl ? mv.val[1] : 0.0

function LinearAlgebra.dot(mv1::Multivector, mv2::Multivector)
    rsbas = intersect(mv1.bas, mv2.bas)
    ln = length(rsbas)
    res = 0.0
    for i in 1:ln
        sn = gaprodsign(rsbas[i], rsbas[i])
        j = findfirst(isequal(rsbas[i]), mv1.bas)
        tp = mv1.val[j] * sn
        j = findfirst(isequal(rsbas[i]), mv2.bas)
        res += tp * mv2.val[j]
    end
    return res
end

#Exponentiation
function Base.exp(a::Multivector)
    s = max(ceil(Int, log(2, dot(a.val, a.val))) - 1, 0)
    a = 1 / 2^s * a
    res = 1 + a
    powa = a
    for i in 2:12
        powa *= a / i
        res += powa
    end
    while s > 0
        res = res * res
        s -= 1
    end
    return res
end

function SimpleGA.bivector_exp(a::Multivector)
    R = exp(project(a, 2))
    delt = R * R' - 1
    return (1 - 0.5 * delt + 0.375 * delt * delt) * R
end

# Approximate comparison
function Base.isapprox(mv1::Multivector, mv2::Multivector)
    tmp = mvtidy(mv1 - mv2)
    return tmp.bas == [basscl] && tmp.val == [0.0]
end
