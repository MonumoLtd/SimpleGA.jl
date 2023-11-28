# Common functions used in testing.

"""Test that off-diagonal basis elements all vanish."""
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

"""Tests that unit vector generators behave correctly over integers"""
function run_test_positive_norm(e1, e2)
    @test e1 * e1 == one(e1 * e2)
    @test e2 * e2 == one(e2 * e1)
    @test 1 + e1 * e2 == e1 * e2 + 1
    @test 1 - e1 * e2 == -e1 * e2 + 1
    @test (1 + e1 * e2) * (1 + e1 * e2) == 2 * e1 * e2
end

function run_test_mixed_norm(e1, f1)
    @test e1 * e1 == one(e1 * f1)
    @test f1 * f1 == -one(f1 * e1)
    @test 1 + e1 * f1 == e1 * f1 + 1
    @test 1 - e1 * f1 == -e1 * f1 + 1
    @test iszero((1 + e1 * f1) * (1 - e1 * f1))
end

"""Test common features of all algebras."""
function run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)
    # Addition
    @test isapprox(1.0 + me1, me1 + 1.0)
    @test isapprox(-1.0 + me1, me1 - 1.0)

    # Distributivity
    @test isapprox(me1 * (me2 + me3), me1 * me2 + me1 * me3)
    @test isapprox(mo1 * (me2 + me3), mo1 * me2 + mo1 * me3)
    @test isapprox(me1 * (mo2 + mo3), me1 * mo2 + me1 * mo3)
    @test isapprox(mo1 * (mo2 + mo3), mo1 * mo2 + mo1 * mo3)

    # Associativity
    @test isapprox(me1 * (me2 * me3), (me1 * me2) * me3)
    @test isapprox(mo1 * (me2 * me3), (mo1 * me2) * me3)
    @test isapprox(me1 * (mo2 * me3), (me1 * mo2) * me3)
    @test isapprox(me1 * (me2 * mo3), (me1 * me2) * mo3)
    @test isapprox(mo1 * (mo2 * me3), (mo1 * mo2) * me3)
    @test isapprox(mo1 * (me2 * mo3), (mo1 * me2) * mo3)
    @test isapprox(me1 * (mo2 * mo3), (me1 * mo2) * mo3)
    @test isapprox(mo1 * (mo2 * mo3), (mo1 * mo2) * mo3)

    # Projection
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

    # Rotation
    R = bivector_exp(v1 * v2)
    no1 = R * mo1 * R'
    no2 = R * mo2 * R'
    @test isapprox(dot(mo1, mo2), dot(no1, no2))
    ne1 = R * me1 * R'
    ne2 = R * me2 * R'
    @test isapprox(dot(me1, me2), dot(ne1, ne2))

    # Reverse
    @test isapprox((me1 + me1') / 2, tr(me1) + project(me1, 4))
    @test isapprox((me1 - me1') / 2, project(me1, 2) + project(me1, 6))
    @test isapprox((mo1 + mo1') / 2, project(mo1, 1) + project(mo1, 5))
    @test isapprox((mo1 - mo1') / 2, project(mo1, 3))
end

"""Tests conversion to Float32"""
function run_conversion_tests(me1, me2, mo1, mo2, E, O)
    @test isapprox(convert(E, me1 * me2), convert(E, me1) * convert(E, me2))
    @test isapprox(convert(E, mo1 * mo2), convert(O, mo1) * convert(O, mo2))
    @test isapprox(convert(O, me1 * mo2), convert(E, me1) * convert(O, mo2))
    @test typeof(convert(E, me1) * convert(E, me2)) == E
    @test typeof(convert(O, mo1) * convert(O, mo2)) == E
    @test typeof(convert(E, me1) * convert(O, mo2)) == O
    cme = convert(E, me1)
    cme32 =
        project(cme, 0) +
        project(cme, 2) +
        project(cme, 4) +
        project(cme, 6) +
        project(cme, 8)
    @test typeof(cme32) == E
    cmo = convert(O, mo1)
    cmo32 = project(cmo, 1) + project(cmo, 3) + project(cmo, 5) + project(cmo, 7)
    @test typeof(cmo32) == O
    @test typeof(cme) == typeof(cme')
    @test typeof(cmo) == typeof(cmo')
end
