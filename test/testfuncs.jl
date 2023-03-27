#Common functions used in testing.

#Function to test off-diagonal basis elements all vanish.
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

#Test shared across all algebras.
function run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)
    #Addition
    @test isapprox(1.0 + me1, me1 + 1.0)
    @test isapprox(-1.0 + me1, me1 - 1.0)

    #Distributivity
    @test isapprox(me1 * (me2 + me3), me1 * me2 + me1 * me3)
    @test isapprox(mo1 * (me2 + me3), mo1 * me2 + mo1 * me3)
    @test isapprox(me1 * (mo2 + mo3), me1 * mo2 + me1 * mo3)
    @test isapprox(mo1 * (mo2 + mo3), mo1 * mo2 + mo1 * mo3)

    #Associativity
    @test isapprox(me1 * (me2 * me3), (me1 * me2) * me3)
    @test isapprox(mo1 * (me2 * me3), (mo1 * me2) * me3)
    @test isapprox(me1 * (mo2 * me3), (me1 * mo2) * me3)
    @test isapprox(me1 * (me2 * mo3), (me1 * me2) * mo3)
    @test isapprox(mo1 * (mo2 * me3), (mo1 * mo2) * me3)
    @test isapprox(mo1 * (me2 * mo3), (mo1 * me2) * mo3)
    @test isapprox(me1 * (mo2 * mo3), (me1 * mo2) * mo3)
    @test isapprox(mo1 * (mo2 * mo3), (mo1 * mo2) * mo3)

    #Projection
    @test isapprox(
        me1,
        project(me1, 0) +
        project(me1, 2) +
        project(me1, 4) +
        project(me1, 6) +
        project(me1, 8),
    )
    @test isapprox(
        mo1, project(mo1, 1) + project(mo1, 3) + project(mo1, 5) + project(mo1, 7)
    )

    #Rotation
    R = bivector_exp(v1 * v2)
    no1 = R * mo1 * R'
    no2 = R * mo2 * R'
    @test isapprox(dot(mo1, mo2), dot(no1, no2))
    ne1 = R * me1 * R'
    ne2 = R * me2 * R'
    @test isapprox(dot(me1, me2), dot(ne1, ne2))

    #Reverse
    @test isapprox((me1 + me1') / 2, tr(me1) + project(me1, 4))
    @test isapprox((me1 - me1') / 2, project(me1, 2) + project(me1, 6))
    @test isapprox((mo1 + mo1') / 2, project(mo1, 1) + project(mo1, 5))
    @test isapprox((mo1 - mo1') / 2, project(mo1, 3))
end
