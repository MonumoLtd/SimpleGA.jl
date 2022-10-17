# Relative tolerance for common comparisons. Loosened slightly from the default due to
# some numerical precision issues.
const ATOL = 1e-14

function test_distributive(
    me1::MultiVector{Even},
    me2::MultiVector{Even},
    me3::MultiVector{Even},
    mo1::MultiVector{Odd},
    mo2::MultiVector{Odd},
    mo3::MultiVector{Odd},
)
    @test me1 * (me2 + me3) ≈ me1 * me2 + me1 * me3 atol=ATOL
    @test mo1 * (me2 + me3) ≈ mo1 * me2 + mo1 * me3 atol=ATOL
    @test me1 * (mo2 + mo3) ≈ me1 * mo2 + me1 * mo3 atol=ATOL
    @test mo1 * (mo2 + mo3) ≈ mo1 * mo2 + mo1 * mo3 atol=ATOL
    return nothing
end

function test_associative(
    me1::MultiVector{Even},
    me2::MultiVector{Even},
    me3::MultiVector{Even},
    mo1::MultiVector{Odd},
    mo2::MultiVector{Odd},
    mo3::MultiVector{Odd},
)
    @test me1 * (me2 * me3) ≈ (me1 * me2) * me3 atol=ATOL
    @test mo1 * (me2 * me3) ≈ (mo1 * me2) * me3 atol=ATOL
    @test me1 * (mo2 * me3) ≈ (me1 * mo2) * me3 atol=ATOL
    @test me1 * (me2 * mo3) ≈ (me1 * me2) * mo3 atol=ATOL
    @test mo1 * (mo2 * me3) ≈ (mo1 * mo2) * me3 atol=ATOL
    @test mo1 * (me2 * mo3) ≈ (mo1 * me2) * mo3 atol=ATOL
    @test me1 * (mo2 * mo3) ≈ (me1 * mo2) * mo3 atol=ATOL
    @test mo1 * (mo2 * mo3) ≈ (mo1 * mo2) * mo3 atol=ATOL
    return nothing
end

function test_rotation(e1::MultiVector{Odd}, e2::MultiVector{Odd}, me::MultiVector{Even})
    f1 = me * e1 * reverse(me)
    f2 = me * e2 * reverse(me)
    @test scalar(f1, f2) ≈ 0 atol=ATOL
    return nothing
end

function test_projection(me::MultiVector{Even})
    @test me ≈ project(me, 0) + project(me, 2) + project(me, 4) + project(me, 6) atol=ATOL
    return nothing
end

function test_projection(mo::MultiVector{Odd})
    @test mo ≈ project(mo, 1) + project(mo, 3) + project(mo, 5) atol=ATOL
    return nothing
end

function test_algebra(e1, e2, me1, me2, me3, mo1, mo2, mo3)
    @testset "distributivity" begin
        test_distributive(me1, me2, me3, mo1, mo2, mo3)
    end

    @testset "associativity" begin
        test_associative(me1, me2, me3, mo1, mo2, mo3)
    end

    @testset "rotation" begin
        test_rotation(e1, e2, me1)
        test_rotation(e1, e2, me2)
        test_rotation(e1, e2, me3)
    end

    @testset "projection" begin
        test_projection(me1)
        test_projection(me2)
        test_projection(me3)
        test_projection(mo1)
        test_projection(mo2)
        test_projection(mo3)
    end

    return nothing
end
