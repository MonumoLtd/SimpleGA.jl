# Common functions used in testing.

"""Return true iff off-diagonal basis elements all vanish."""
function _do_off_diagonal_elements_vanish(bas)
    n = length(bas)
    for i in 1:(n - 1)
        for j in (i + 1):n
            ok = isequal(bas[i] * bas[j] + bas[j] * bas[i], zero(bas[i] * bas[j]))
            ok || return false
        end
    end
    return true
end

"""Run common tests given a basis."""
function run_basis_tests(basis)
    @test _do_off_diagonal_elements_vanish(basis)

    @testset "promotion" begin
        # If our elements are not bitstypes, then we cannot compare values with identity. 
        # This will be the case for GA(4, 4) and GA(32, 32).
        # For these algebras, we will just compare the _types_ with identity when checking
        # for promotion.
        # NOTE: In practice it is probably good enough to always just check the types, since
        #   we later perform an equality check on values.
        can_test_identity = all.(isbits(basis))

        # Tests for https://github.com/MonumoLtd/SimpleGA.jl/issues/18
        for e in basis
            @test 1.0 * e !== e
            @test 1.0f0 * e !== 1.0 * e
            if can_test_identity
                @test promote(1.0f0 * e, e) === (1.0f0 * e, 1.0f0 * e)
                @test promote(1.0 * e, e) === (1.0 * e, 1.0 * e)
                @test promote(1.0 * e, 1.0f0 * e) === (1.0 * e, 1.0 * e)
            else
                @test typeof(promote(1.0f0 * e, e)) === typeof((1.0f0 * e, 1.0f0 * e))
                @test typeof(promote(1.0 * e, e)) === typeof((1.0 * e, 1.0 * e))
                @test typeof(promote(1.0 * e, 1.0f0 * e)) === typeof((1.0 * e, 1.0 * e))
            end

            # If promotion works, so should equality testing
            @test 1.0 * e == e
            @test isequal(1.0 * e, e)
            @test isapprox(1.0 * e, e)
        end

        # Check that we also promote even elements as expected.
        @test length(basis) >= 2

        # Construct a bivector B; this will be even.
        B = basis[1] * basis[2]
        @test project(B, 2) == B  # Convince ourselves that this is indeed a bivector
        @test 1.0 * B !== B
        @test 1.0f0 * B !== 1.0 * B
        if can_test_identity
            @test promote(1.0f0 * B, B) === (1.0f0 * B, 1.0f0 * B)
            @test promote(1.0 * B, B) === (1.0 * B, 1.0 * B)
            @test promote(1.0 * B, 1.0f0 * B) === (1.0 * B, 1.0 * B)
        else
            @test typeof(promote(1.0f0 * B, B)) === typeof((1.0f0 * B, 1.0f0 * B))
            @test typeof(promote(1.0 * B, B)) === typeof((1.0 * B, 1.0 * B))
            @test typeof(promote(1.0 * B, 1.0f0 * B)) === typeof((1.0 * B, 1.0 * B))
        end
    end
end

"""Run common tests based on a multivector type.

Note that Even and Odd _might_ be the same type.
"""
function test_type(Even::Type, Odd::Type)
    # Scalar construction
    ae = zero(Even{Int8})
    @test zero(ae) == ae
    ao = zero(Odd{Int8})
    @test zero(ao) == ao

    be = one(Even{Int8})
    # NB: There is no concept of "one" for Odd multivectors

    # arrays of ones and zeros should work too.
    aes = zeros(Even{Int8}, 3, 3)
    @test all(aes .== ae)
    aos = zeros(Odd{Int8}, 3, 3)
    @test all(aos .== ao)

    bes = ones(Even{Int8}, 3, 3)
    @test all(bes .== be)
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
    # isapprox kwargs. Tests for
    #   https://github.com/MonumoLtd/SimpleGA.jl/issues/20
    @test isapprox(me1, me1)
    @test isapprox(me1, me1; rtol=1e-5)
    @test !isapprox(me1, me2)
    @test !isapprox(me1, me2; rtol=1e-5)

    # Test that we can compare arrays of multivectors
    @test isapprox([me1], [me1])
    @test isapprox([me1], [me1]; rtol=1e-5)
    @test !isapprox([me1], [me2])
    @test !isapprox([me1], [me2]; rtol=1e-5)
    @test isapprox([me1;;], [me1;;])
    @test isapprox([me1;;], [me1;;]; rtol=1e-5)
    @test !isapprox([me1;;], [me2;;])
    @test !isapprox([me1;;], [me2;;]; rtol=1e-5)

    # Test that isapprox is false for comparisons between odd and even elements
    @test !isapprox(me1, mo1)
    @test !isapprox(me1, mo1; rtol=1e-5)
    @test !isapprox([me1], [mo1])
    @test !isapprox([me1;;], [mo1;;])
    @test !isapprox([me1], [mo1]; rtol=1e-5)
    @test !isapprox([me1;;], [mo1;;]; rtol=1e-5)

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
