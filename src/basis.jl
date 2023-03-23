"""
The core interface function Basis.
"""

function basis(p::Integer)
    return if p <= 2
        bas20
    elseif p == 3
        bas30
    elseif p == 4
        bas40
    else
        bas64
    end
end

function basis(p::Integer, q::Integer)
    if q == 0
        return basis(p)
    elseif q == 1
        p <= 4 ? basCGA : bas64
    elseif q == 2
        if p <= 1
            basSTA
        elseif p <= 3
            bas33
        elseif p == 4
            bas44
        else
            bas64
        end
    elseif q == 3
        if p <= 1
            basSTA
        elseif p == 2
            bas24
        elseif p == 3
            bas33
        elseif p == 4
            bas44
        else
            bas64
        end
    elseif q == 4
        p <= 4 ? bas44 : bas64
    else
        return bas64
    end
end

function basis(p::Integer, q::Integer, r::Integer)
    if r == 0
        return basis(p, q)
    elseif r == 1
        p <= 3 && q == 0 ? basPGA : basis(p + 1, q + 1)
    else
        return basis(p + r, q + r)
    end
end

#Included for convenience for people who prefer to refer to an algebra by name.
function basis(alg::String)
    if alg == "GA20"
        return bas20
    elseif alg == "GA30"
        return bas30
    elseif alg == "GA40"
        return bas40
    elseif alg == "PGA"
        return basPGA
    elseif alg == "CGA"
        return basCGA
    elseif alg == "STA"
        return basSTA
    elseif alg == "GA33"
        return bas33
    elseif alg == "GA24"
        return bas24
    elseif alg == "GA44"
        return bas44
    elseif alg == "GA64"
        return bas64
    else throw(ArgumentError(alg * " is not a recognised algebra."))
    end
end


#Test functions

function testbas(bas)
    res = true
    n = length(bas)
    for i in 1:(n - 1)
        for j in (i + 1):n
            res = isequal(res && bas[i] * bas[j] + bas[j] * bas[i], zero(bas[i] * bas[j]))
        end
    end
    return res
end
