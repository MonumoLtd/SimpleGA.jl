# Short-cuts for the algebra we wish to use.
e1 = GA20.e1
e2 = GA20.e2
I2 = GA20.I2

@testset "stringiness" begin
    scalar = GA20.MV{Even}(42.0)
    @test string(scalar) == "42.0"
    @test string(e1) == "1.0 e1"
    @test string(e2) == "1.0 e2"
    @test string(I2) == "1.0 I2"

    @test string(-e1) == "-1.0 e1"
    @test string(-e2) == "-1.0 e2"
    @test string(-I2) == "-1.0 I2"

    @test string(42 + 3 * I2) == "42.0 + 3.0 I2"
    @test string(42 + I2 * 3) == "42.0 + 3.0 I2"
end

@testset "no parity mixing" begin
    @test_throws MethodError 42 + 3 * e1
    @test_throws MethodError 42 + 3 * e2
    @test_throws MethodError e2 + I2 * 1.5
end

# Run common tests on the algebra
me1 = rand() + rand() * e1 * e2
me2 = rand() + rand() * e1 * e2
me3 = rand() + rand() * e1 * e2
mo1 = rand() * e1 + rand() * e2
mo2 = rand() * e1 + rand() * e2
mo3 = rand() * e1 + rand() * e2
test_algebra(e1, e2, me1, me2, me3, mo1, mo2, mo3)
