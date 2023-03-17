#= 
Code for GA(32,32).
Not optimised in any way, but useful for some experimentation.
=#

module GA64

using LinearAlgebra

include("GAcore64.jl")
import Base.show



bldbas = [
0x0000000000000001,0x0000000000000002,0x0000000000000007,0x000000000000000b,
0x000000000000001f,0x000000000000002f,0x000000000000007f,0x00000000000000bf,
0x00000000000001ff,0x00000000000002ff,0x00000000000007ff,0x0000000000000bff,
0x0000000000001fff,0x0000000000002fff,0x0000000000007fff,0x000000000000bfff,
0x000000000001ffff,0x000000000002ffff,0x000000000007ffff,0x00000000000bffff,
0x00000000001fffff,0x00000000002fffff,0x00000000007fffff,0x0000000000bfffff,
0x0000000001ffffff,0x0000000002ffffff,0x0000000007ffffff,0x000000000bffffff,
0x000000001fffffff,0x000000002fffffff,0x000000007fffffff,0x00000000bfffffff,
0x00000001ffffffff,0x00000002ffffffff,0x00000007ffffffff,0x0000000bffffffff,
0x0000001fffffffff,0x0000002fffffffff,0x0000007fffffffff,0x000000bfffffffff,
0x000001ffffffffff,0x000002ffffffffff,0x000007ffffffffff,0x00000bffffffffff,
0x00001fffffffffff,0x00002fffffffffff,0x00007fffffffffff,0x0000bfffffffffff,
0x0001ffffffffffff,0x0002ffffffffffff,0x0007ffffffffffff,0x000bffffffffffff,
0x001fffffffffffff,0x002fffffffffffff,0x007fffffffffffff,0x00bfffffffffffff,
0x01ffffffffffffff,0x02ffffffffffffff,0x07ffffffffffffff,0x0bffffffffffffff,
0x1fffffffffffffff,0x2fffffffffffffff,0x7fffffffffffffff,0xbfffffffffffffff
]

bas64 = map(n->Multivector([n],[1.0]), bldbas)
export bas64, construct64

struct Blade
    bas::UInt64
    val::Number
end



function Base.isless(x::Blade,y::Blade)
    if grd(x.bas) < grd(y.bas)
        return true
    elseif grd(x.bas) > grd(y.bas)
        return false
    else
        return isless(x.bas,y.bas)
    end
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



function vece(n::Int)
    if n == 0
        return Multivector([0x0000000000000000],[1.0])
    end
    #Keep in range 1..32
    nn = ((((n-1) % 32) + 32) % 32) +1 
    return Multivector([bldlut[2*n-1]],[1.0])
end


function vece(n::Int, val::Float64)
    return val*vece(n)
end

function vecf(n::Int)
    if n == 0
        return Multivector([0x0000000000000000],[1.0])
    end
    #Keep in range 1..32
    nn = ((((n-1) % 32) + 32) % 32) +1 
    return Multivector([bldlut[2*n]],[1.0])
end

function vecf(n::Int, val::Float64)
    return val*vecf(n)
end

function bdptype(nn::UInt64)
    bldn = bldconvert(nn)
    res=""
    for i in 1:32 
        if trailing_zeros(bldn) == 0
            res = res *"e"* string(i)
        end
        bldn = bldn >> 1
        if trailing_zeros(bldn) == 0
            res = res *"f"* string(i)
        end
        bldn = bldn >> 1
    end
    return res
end

function mvtype(mv::Multivector)
    blds = mvtoblds(mv)
    sort!(blds)
    res = string(blds[1].val) * bdptype(blds[1].bas)
    if length(blds) == 1
        return res
    end
    for i in 2:length(blds)
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