#=
Core code for the implementation of GA(4,4).
Uses a Uint8 representation of basis blades and bitwise operations.
Used for checking other algebras.
=#

#This avoids rounding error bloating multivectors. Assumes using FP64. Set to zero if unsure
const mvtol = 1e-14
sparsify(x, eps) = abs(x) < eps ? 0.0 : x

#The Multivector type assumes that the blade list is unique and in order. But we want to avoid checking this at runtime.
#Only use this constructor if you are certain the blade list is correct. If not, use construct44()

struct Multivector{T<:Real} <: Number
    bas::Vector{UInt8}
    val::Vector{T}
end

Multivector(ns, vs) = Multivector{typeof(vs[1])}(convert(Vector{UInt8}, ns), vs)

# NOTE: we do not define equality in terms of the individual components, since we have many ways to represent "zero" that
#   we wish to be equivalent.
Base.:(==)(mv1::Multivector, mv2::Multivector) = iszero((mv1 - mv2).val)
Base.hash(a::Multivector, h::UInt) = hash(a.bas, hash(a.val, hash(:Multivector, h)))

function construct44(T, bs, vs)
    if length(bs) != length(unique(bs))
        error("List of blades must be unique")
    else
        p = sortperm(bs)
        return Multivector{T}(sort(bs), vs[p])
    end
end

function Base.convert(::Type{Multivector{T}}, a::Multivector) where {T<:Real}
    return Multivector{T}(a.bas, convert.(T, a.val))
end
function Base.promote_rule(
    ::Type{Multivector{S}}, ::Type{Multivector{T}}
) where {S<:Real,T<:Real}
    return Multivector{promote_type(S, T)}
end
Base.zero(::Type{Multivector{T}}) where {T} = Multivector([0x00], [zero(T)])
Base.one(::Type{Multivector{T}}) where {T} = Multivector([0x00], [one(T)])
Base.zero(a::Multivector) = zero(typeof(a))
Base.one(a::Multivector) = one(typeof(a))

#Addition / subtraction
Base.:(-)(mv::Multivector) = Multivector(mv.bas, -mv.val)

function Base.:(+)(mv1::Multivector, mv2::Multivector)
    res =
        SparseVector(256, mv1.bas .+ 1, mv1.val) + SparseVector(256, mv2.bas .+ 1, mv2.val)
    sps = sparse(sparsify.(res, mvtol))
    return length(sps.nzind) == 0 ? zero(mv1) : Multivector(sps.nzind .- 1, sps.nzval)
end

Base.:(+)(nm::Real, mv::Multivector) = mv + Multivector([0x00], [nm])
Base.:(+)(mv::Multivector, nm::Real) = nm + mv
Base.:(-)(nm::Real, mv::Multivector) = nm + (-mv)
Base.:(-)(mv::Multivector, nm::Real) = (-nm) + mv
Base.:(-)(mv1::Multivector, mv2::Multivector) = mv1 + (-(mv2))

#Multiplication
function gaprodsign(bld1::UInt8, bld2::UInt8)
    tp1 = xor(bld2, bld2 << 1)
    cntones = count_ones((bld1 & 0xaa) & tp1)
    return convert(Int8, 1 - 2 * (cntones % 2))
end

function Base.:(*)(mv1::Multivector, mv2::Multivector)
    res = zeros(typeof(mv1.val[1]), 256)
    for i in 1:length(mv1.bas)
        for j in 1:length(mv2.bas)
            res[xor(mv1.bas[i], mv2.bas[j]) + 1] +=
                gaprodsign(mv1.bas[i], mv2.bas[j]) * mv1.val[i] * mv2.val[j]
        end
    end
    sps = sparse(sparsify.(res, mvtol))
    return length(sps.nzind) == 0 ? zero(mv1) : Multivector(sps.nzind .- 1, sps.nzval)
end

Base.:(*)(num::Real, mv::Multivector) = Multivector(mv.bas, num * mv.val)
Base.:(*)(mv::Multivector, num::Real) = num * mv
Base.:(/)(mv::Multivector, num::Real) = (1 / num) * mv

#Reverse

function grd(xin::UInt8)
    if count_ones(xin & 0xc0) == 1
        xin = xor(xin, 0x3f)
    end
    if count_ones(xin & 0x30) == 1
        xin = xor(xin, 0x0f)
    end
    if count_ones(xin & 0x0c) == 1
        xin = xor(xin, 0x03)
    end
    return count_ones(xin)
end

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
    if iszero(ln)
        return zero(mv)
    end
    rsval = zeros(typeof(mv.val[1]), ln)
    for i in 1:ln
        j = findfirst(isequal(rsbas[i]), mv.bas)
        rsval[i] = mv.val[j]
    end
    return Multivector(rsbas, rsval)
end

LinearAlgebra.tr(mv::Multivector) = iszero(mv.bas[1]) ? mv.val[1] : 0.0

function LinearAlgebra.dot(mv1::Multivector, mv2::Multivector)
    rsbas = intersect(mv1.bas, mv2.bas)
    ln = length(rsbas)
    res = zero(mv1.val[1])
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

# XXX: Presumably the default for `atol` here is required due to how we are
#   representing the multivector?
function Base.isapprox(mv1::Multivector, mv2::Multivector; atol=1e-6, kwargs...)
    return isapprox(
        SparseVector(256, mv1.bas .+ 1, mv1.val),
        SparseVector(256, mv2.bas .+ 1, mv2.val);
        atol,
        kwargs...,
    )
end
