#= 
Implementation of GA(32,32) in Julia using Uint64 bitwise operations. 
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


#The Multivector type assumes that the blade list is unique and in order. But we want to avoid checking this at runtime.
#Only use this constructor if you are certain the blade list is correct. If not, use construct64()
struct Multivector
    bas::Vector{UInt64}
    val::Vector{Float64}
end

basscl = 0x0000000000000000
mvzero = Multivector([basscl], [0.0])
mvone = Multivector([basscl], [1.0])


function construct64(bs,vs)
    if length(bs) != length(unique(bs))
        error("List of blades must be unique")
    else
        p = sortperm(bs)
        return Multivector(sort(bs),vs[p])
    end
end


#This avoids rounding error bloating multivectors. 
const mvtol = 1e-14

function mvtidy(mv::Multivector)
    ln=length(filter(x->!isapprox(x,0.0;atol = mvtol), mv.val))
    if ln==0
        return mvzero
    end
    rsbas = zeros(UInt64,ln)
    rsval = zeros(Float64,ln)
    j=1
    for i = 1:length(mv.bas)
        if !isapprox(mv.val[i],0.0)
            rsbas[j]=mv.bas[i]
            rsval[j]=mv.val[i]
            j+=1
        end
    end
    return Multivector(rsbas,rsval)
end


#Addition / subtraction
function -(mv::Multivector)
    Multivector(mv.bas,-mv.val)
end


function +(mv1::Multivector, mv2::Multivector)
    rsbas=sort(union(mv1.bas,mv2.bas))
    l1 = length(mv1.bas)
    l2 = length(mv2.bas)
    ln=length(rsbas)
    rsval=zeros(Float64,ln)
    i = 1
    j = 1
    for k in 1:ln
        if i <= l1 && rsbas[k] == mv1.bas[i]
            rsval[k]+=mv1.val[i]
            i+=1
        end
        if j <= l2 && rsbas[k] == mv2.bas[j]
            rsval[k]+=mv2.val[j]
            j+=1
        end
    end
    return Multivector(rsbas, rsval)
end


function +(nm::Number,mv::Multivector)
    mv2 = Multivector([basscl],[convert(Float64,nm)])
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

function gaprodsign(bld1, bld2)
    tp1 = xor( bld2, bld2 << 1 )
    cntones = count_ones((bld1 & 0xaaaaaaaaaaaaaaaa) & tp1)
    return convert(Int8, 1 - 2* (cntones % 2))
end



function *(num::Number,mv::Multivector)
    num = convert(Float64,num)
    return Multivector(mv.bas,num*mv.val)
end


function *(mv1::Multivector,mv2::Multivector)
    res = 0
    l1 = length(mv1.bas)
    for i in 1:l1
        bs = map(bd->xor(mv1.bas[i],bd), mv2.bas)
        vs = map((bd,vl) -> mv1.val[i]*vl*gaprodsign(mv1.bas[i],bd),mv2.bas,mv2.val)
        p = sortperm(bs)
        res += Multivector(sort(bs),vs[p])
    end
    return mvtidy(res)
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

#Converts conformal split to traditional rep for grade operation
function bldconvert(xin::UInt64) 
    res = xin
    msk = 0xc000000000000000
    chk = ~msk
    for i in 1:31
        if count_ones(res & msk) == 1
            res = xor(res,chk)
        end
        msk = msk >> 2
        chk = chk >> 2
    end
    return res
end

function grd(xin::UInt64)
    return count_ones(bldconvert(xin))
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
    if mv.bas[1]==basscl
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

#Comparison
function Base.isapprox(mv1::Multivector, mv2::Multivector)
    tmp = mvtidy(mv1-mv2)
    return tmp.bas == [basscl] && tmp.val ==[0.0]
end