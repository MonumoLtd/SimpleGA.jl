# Short-cuts for the algebra we wish to use.
e1 = GA30.e1
e2 = GA30.e2
e3 = GA30.e3
I3 = GA30.I3

@testset "stringiness" begin
    # FIXME need a nicer way of creating a scalar across different algebras.
    scalar = GA30.MV{Even}(42.0, 0.0, 0.0, 0.0)
    @test string(scalar) == "42.0"
    @test string(e1) == "1.0 e1"
    @test string(e2) == "1.0 e2"
    @test string(e3) == "1.0 e3"
    @test string(I3) == "1.0 I3"

    @test string(-e1) == "-1.0 e1"
    @test string(-e2) == "-1.0 e2"
    @test string(-e3) == "-1.0 e3"
    @test string(-I3) == "-1.0 I3"

    @test string(42 + 3 * e1 * e2) == "42.0 + 3.0 e1e2"
    @test string(42 + e1 * 3 * e2) == "42.0 + 3.0 e1e2"
    @test string(42 * e1 + 5 * I3) == "42.0 e1 + 5.0 I3"
end

@testset "no parity mixing" begin
    @test_throws MethodError 42 + 3 * e1
    @test_throws MethodError 42 + 3 * e2
    @test_throws MethodError e2 + I2 * 1.5
end

# Run common tests on the algebra
me1 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
me2 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
me3 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
mo1 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1
mo2 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1
mo3 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1

test_algebra(e1, e2, me1, me2, me3, mo1, mo2, mo3)
