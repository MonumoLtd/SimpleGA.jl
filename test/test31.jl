#Test suite for GA31.
#Test stand-alone results and compares with GA(4,4)

basGA31 = GA31.basis
(e1, e2, e3, f3) = (basGA31[1], basGA31[2], basGA31[3], basGA31[4])
I4 = e1 * e2 * e3 * f3

@test map(x -> dot(x, x), basGA31) == [1, 1, 1, -1]
@test testbas(basGA31)

#! format:off
me1 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 + f3*(rand()*e1 - rand()*e2 +rand()*e3) + I4*rand()
me2 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 + f3*(rand()*e1 - rand()*e2 +rand()*e3) + I4*rand()
me3 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 + f3*(rand()*e1 - rand()*e2 +rand()*e3) + I4*rand()
mo1 = rand()*e1 + rand()*e2 + e3*rand() + f3/rand() + I4*(rand()*e1 + e2*rand() + e3*rand() - rand()*f3)
mo2 = rand()*e1 + rand()*e2 + e3*rand() + f3/rand() + I4*(rand()*e1 + e2*rand() + e3*rand() - rand()*f3)
mo3 = rand()*e1 + rand()*e2 + e3*rand() + f3/rand() + I4*(rand()*e1 + e2*rand() + e3*rand() - rand()*f3)

#! format:on

#Comparison with GA(4,4)
bas44 = GA44.basis
(E1, E2, E3, F3) = (bas44[1], bas44[2], bas44[3], bas44[8])
arr1 = rand(4)
v1 = inject(arr1, basGA31)
V1 = inject(arr1, [E1, E2, E3, F3])
arr2 = rand(4)
v2 = inject(arr2, basGA31)
V2 = inject(arr2, [E1, E2, E3, F3])
arr3 = rand(4)
v3 = inject(arr3, basGA31)
V3 = inject(arr3, [E1, E2, E3, F3])
arr4 = rand(4)
v4 = inject(arr4, basGA31)
V4 = inject(arr4, [E1, E2, E3, F3])
arr5 = rand(4)
v5 = inject(arr5, basGA31)
V5 = inject(arr5, [E1, E2, E3, F3])
@test isapprox(dot(v1, v1), dot(V1, V1))
@test isapprox(dot(v1 * v2 * v3, e1), dot(V1 * V2 * V3, E1))
@test isapprox(dot(v1 * v2 * v3, e2), dot(V1 * V2 * V3, E2))
@test isapprox(dot(v1 * v2 * v3, e3), dot(V1 * V2 * V3, E3))
@test isapprox(dot(v1 * v2 * v3, f3), dot(V1 * V2 * V3, F3))
@test isapprox(embed(v1 * v2 * v3), V1 * V2 * V3)
@test isapprox(embed(v1 * v2 * v3 * v4), V1 * V2 * V3 * V4)
@test isapprox(embed(v1 * v2 * v3 * v4 * v5), V1 * V2 * V3 * V4 * V5)
@test isapprox(embed(exp(v1 * v2)), exp(V1 * V2))
@test isapprox(embed(bivector_exp(v1 * v2)), bivector_exp(V1 * V2))

run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)

# Conversion
run_conversion_tests(me1, me2, mo1, mo2, GA31.Even{Float32}, GA31.Odd{Float32})