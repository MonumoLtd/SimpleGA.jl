"""
The core interface function Basis.
"""

function basis(p::Integer)
    return if p <= 2
        GA20.basis
    elseif p == 3
        GA30.basis
    elseif p == 4
        GA40.basis
    else
        GA3232.basis
    end
end

function basis(p::Integer, q::Integer)
    return if q == 0
        basis(p)
    elseif q == 1
        p <= 4 ? CGA.basis : GA3232.basis
    elseif q == 2
        if p <= 1
            STA.basis
        elseif p <= 3
            GA33.basis
        elseif p == 4
            GA44.basis
        else
            GA3232.basis
        end
    elseif q == 3
        if p <= 1
            STA.basis
        elseif p == 2
            GA24.basis
        elseif p == 3
            GA33.basis
        elseif p == 4
            GA44.basis
        else
            GA3232.basis
        end
    elseif q == 4
        if p <= 2
            GA24.basis
        else
            p <= 4 ? GA44.basis : GA3232.basis
        end
    else
        GA3232.basis
    end
end

function basis(p::Integer, q::Integer, r::Integer)
    return if r == 0
        basis(p, q)
    elseif r == 1
        p <= 3 && q == 0 ? PGA.basis : basis(p + 1, q + 1)
    else
        basis(p + r, q + r)
    end
end

function basis()
    println("Supported algebras and names ")
    println("GA(2,0); GA20")
    println("GA(3,0); GA30")
    println("GA(4,0); GA40")
    println("GA(1,3); STA")
    println("GA(3,0,1); PGA")
    println("GA(4,1); CGA")
    println("GA(2,4); GA24")
    println("GA(3,3); GA33")
    println("GA(4,4); GA44")
    println("GA(32,32); GA3232")
end
