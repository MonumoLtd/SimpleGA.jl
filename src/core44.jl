#= 
Core code for the implementation of GA(4,4).
Uses a Uint8 representation of basis blades and bitwise operations. 
Used for checking other algebras.
=#

import Base.:*
import Base.:+
import Base.:-
import Base.:/
import Base.exp
import LinearAlgebra.tr
import LinearAlgebra.dot
import LinearAlgebra.adjoint
import ..project
import ..expb

using SparseArrays

#This avoids rounding error bloating multivectors. Set to zero if unsure
const mvtol = 1e-14
sparsify(x, eps) = abs(x) < eps ? 0.0 : x

#The Multivector type assumes that the blade list is unique and in order. But we want to avoid checking this at runtime.
#Only use this constructor if you are certain the blade list is correct. If not, use construct44()
struct Multivector
    bas::Vector{UInt8}
    val::Vector{Float64}
end

function construct44(bs,vs)
    if length(bs) != length(unique(bs))
        error("List of blades must be unique")
    else
        p = sortperm(bs)
        return Multivector(sort(bs),vs[p])
    end
end

mvzero = Multivector([0x00], [0.0])

#Addition / subtraction
function -(mv::Multivector)
    Multivector(mv.bas,-mv.val)
end


function +(mv1::Multivector, mv2::Multivector)
    res = SparseVector(256,mv1.bas .+ 1, mv1.val) + SparseVector(256,mv2.bas .+ 1, mv2.val)
    sps = sparse(sparsify.(res,mvtol))
    return Multivector(sps.nzind .- 1,sps.nzval)
end


function +(nm::Number,mv::Multivector)
    mv2 = Multivector([0x00],[convert(Float64,nm)])
    return (mv+mv2)
end

function +(mv::Multivector,nm::Number,)
    return nm+mv
end

function -(nm::Number,mv::Multivector)
    return nm+(-mv)
end

function -(mv::Multivector,nm::Number,)
    return (-nm)+mv
end

function -(mv1::Multivector,mv2::Multivector)
    return mv1+(-(mv2))
end

#Multiplication
function gaprodsign(bld1::UInt8, bld2::UInt8)
    #Expects UInt8s
    tp1 = xor( bld2, bld2 << 1 )
    cntones = count_ones((bld1 & 0xaa) & tp1)
    return convert(Int8, 1 - 2* (cntones % 2))
end


function *(mv1::Multivector,mv2::Multivector)
    res = zeros(Float64,256)
    for i in 1:length(mv1.bas)
        for j in 1:length(mv2.bas)
            res[xor(mv1.bas[i],mv2.bas[j])+1] += gaprodsign(mv1.bas[i],mv2.bas[j])*mv1.val[i]*mv2.val[j]
        end
    end
    sps = sparse(sparsify.(res,mvtol))
    return Multivector(sps.nzind .- 1,sps.nzval)
end


function *(num::Number,mv::Multivector)
    num = convert(Float64,num)
    return Multivector(mv.bas,num*mv.val)
end

function *(mv::Multivector,num::Number)
    return num*mv
end

#Division by a real
function /(mv::Multivector,num::Number)
    num = convert(Float64,1/num)
    return num*mv
end


#Reverse

function grd(xin::UInt8)
    if count_ones(xin & 0xc0) == 1
        xin = xor(xin,0x3f)
    end
    if count_ones(xin & 0x30) == 1
        xin = xor(xin,0x0f)
    end
    if count_ones(xin & 0x0c) == 1
        xin = xor(xin, 0x03)
    end
    return count_ones(xin)
end

function adjoint(mv::Multivector)
    rsval = similar(mv.val)
    for i in 1:length(mv.bas)
        if isodd(div(grd(mv.bas[i]),2))
            rsval[i]=-mv.val[i]
        else
            rsval[i]= mv.val[i]
        end
    end
    return Multivector(mv.bas,rsval)
end


#Grade and projection
function project(mv::Multivector,n::Int64)
    rsbas = filter(x->grd(x)==n,mv.bas)
    ln = length(rsbas)
    if ln == 0
        return mvzero
    end
    rsval = zeros(Float64,ln)
    for i in 1:ln
        j=findfirst(isequal(rsbas[i]),mv.bas)
        rsval[i]=mv.val[j]
    end
    return Multivector(rsbas, rsval)
end

function tr(mv::Multivector)
    if mv.bas[1]==0x00
        return mv.val[1]
    else
        return 0.0
    end
end

function dot(mv1::Multivector,mv2::Multivector)
    rsbas=intersect(mv1.bas,mv2.bas)
    ln=length(rsbas)
    res=0.0
    for i in 1:ln
        sn=gaprodsign(rsbas[i], rsbas[i])
        j=findfirst(isequal(rsbas[i]),mv1.bas)
        tp = mv1.val[j]*sn
        j=findfirst(isequal(rsbas[i]),mv2.bas)
        res += tp*mv2.val[j]
    end
    return res
end


#Exponentiation
function exp(a::Multivector)
    s = max(ceil(Int,log(2,dot(a.val,a.val)))-1,0)
    a = 1/2^s*a
    res = 1+a
    powa = a
    for i in 2:12
        powa *= a/i
        res += powa
    end
    while s > 0
        res = res*res
        s -= 1
    end
    return res
end

function expb(a::Multivector)
    R = exp(project(a,2))
    delt = R*R'-1
    return (1-0.5*delt + 0.375*delt*delt)*R
end


Base.isapprox(mv1::Multivector, mv2::Multivector) = isapprox(SparseVector(256,mv1.bas .+ 1, mv1.val), 
        SparseVector(256,mv2.bas .+ 1, mv2.val))

