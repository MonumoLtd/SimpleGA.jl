#= 
Code for GA(4,4)
=#


module GA44

using LinearAlgebra

include("GA44Core.jl")
import Base.show


struct Blade
    bas::UInt8
    val::Number
end


function bldless(x::Blade,y::Blade)
    if grd(x.bas) < grd(y.bas)
        return true
    elseif grd(x.bas) > grd(y.bas)
        return false
    else
        return isless(x.bas,y.bas)
    end
end


# Constructors
const e1 = Multivector([parse(UInt8, "00000001", base=2)],[1.0])
const e2 = Multivector([parse(UInt8, "00000111", base=2)],[1.0])
const e3 = Multivector([parse(UInt8, "00011111", base=2)],[1.0])
const e4 = Multivector([parse(UInt8, "01111111", base=2)],[1.0])
const f1 = Multivector([parse(UInt8, "00000010", base=2)],[1.0])
const f2 = Multivector([parse(UInt8, "00001011", base=2)],[1.0])
const f3 = Multivector([parse(UInt8, "00101111", base=2)],[1.0])
const f4 = Multivector([parse(UInt8, "10111111", base=2)],[1.0])

bas44 = [e1,e2,e3,e4,f1,f2,f3,f4]
export bas44, construct44

#Additional functions

#Removes zeros from a multivector.
#Used for pretty typing. Can be used (with care) to optimise.
function mvtidy(mv::Multivector)
    ln=length(filter(x->!isapprox(x,0.0;atol = 1e-12), mv.val))
    if ln==0
        return Multivector([0x00],[0.0])
    end
    rsbas = zeros(UInt8,ln)
    rsval = zeros(Float64,ln)
    j=1
    for i = 1:length(mv.bas)
        if !isapprox(mv.val[i],0.0;atol = 1e-12)
            rsbas[j]=mv.bas[i]
            rsval[j]=mv.val[i]
            j+=1
        end
    end
    return Multivector(rsbas,rsval)
end


function bdptype(nn::UInt8)
    res=""
    tp = isodd(count_ones( nn & 0b11111101)) ? "e1" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b11111110)) ? "f1" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b11110100)) ? "e2" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b11111000)) ? "f2" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b11010000)) ? "e3" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b11100000)) ? "f3" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b01000000)) ? "e4" : ""
    res *= tp
    tp = isodd(count_ones( nn & 0b10000000)) ? "f4" : ""
    res *= tp
    return res
end


function mvtoblds(mvin::Multivector)
    mv=mvtidy(mvin)
    ln = length(mv.bas)
    res=Array{Blade}(undef,ln)
    for i in 1:ln
        res[i]=Blade(mv.bas[i],mv.val[i])
    end
    return res
end

function mvtype(mv::Multivector)
    blds = mvtoblds(mv)
    sort!(blds,lt=bldless)
    res = string(blds[1].val) * bdptype(blds[1].bas)
    n = length(blds)
    if n == 1
        return res
    end
    for i in 2:n
        tp = " + " * string(blds[i].val) * bdptype(blds[i].bas)
        res *= tp
    end
    return res
end

function Base.show(io::IO, mv::Multivector)
    print(mvtype(mv))
end


function Base.show(io::IO, ::MIME"text/plain", mvs::Vector{Multivector})
    n= length(mvs)
    println(io,n,"-element Vector{Multivector}")
    for i in eachindex(mvs)
    println(io, " ", mvtype(mvs[i]))
    end
end



end #Module