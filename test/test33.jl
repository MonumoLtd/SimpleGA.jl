#Test suite for GA(3,3).
#Test stand-alone results and compares with GA(4,4)

bas33 = GA33.basis
(e1, e2, e3, f1, f2, f3) = (bas33[1], bas33[2], bas33[3], bas33[4], bas33[5], bas33[6])

@test map(x -> dot(x, x), bas33) == [1, 1, 1, -1, -1, -1]
@test testbas(bas33)

#! format:off
me1 = rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2 +
    e3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*f3*(rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2)
me2 = rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2 +
    e3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*f3*(rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2)
me3 = rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2 +
    e3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*f3*(rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2)
mo1 = (rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    f3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    e3*f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ))
mo2 = (rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    f3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    e3*f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ))
mo3 = (rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ) ) +
    e3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    f3*( rand() + rand()*e1*e2 + e1*f1*rand() + e1*rand()*f2 - rand()*f1*e2 + e2*f2/rand() + rand()*f1*f2 + rand()*e1*e2*f1*f2) +
    e3*f3*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() + e1*e2*f1*f2*(rand()*e1 + rand()*e2 + f1*rand() + f2*rand() ))
#! format:on

#Comparison with GA(4,4)
bas44 = GA44.basis
(E1, E2, E3, F1, F2, F3) = (bas44[1], bas44[2], bas44[3], bas44[5], bas44[6], bas44[7])
arr1 = rand(6)
v1 = inject(arr1, bas33)
V1 = inject(arr1, [E1, E2, E3, F1, F2, F3])
arr2 = rand(6)
v2 = inject(arr2, bas33)
V2 = inject(arr2, [E1, E2, E3, F1, F2, F3])
arr3 = rand(6)
v3 = inject(arr3, bas33)
V3 = inject(arr3, [E1, E2, E3, F1, F2, F3])
arr4 = rand(6)
v4 = inject(arr4, bas33)
V4 = inject(arr4, [E1, E2, E3, F1, F2, F3])
@test isapprox(dot(v1, v1), dot(V1, V1))
@test isapprox(dot(v1 * v2 * v3, e1), dot(V1 * V2 * V3, E1))
@test isapprox(dot(v1 * v2 * v3, e2), dot(V1 * V2 * V3, E2))
@test isapprox(dot(v1 * v2 * v3, e3), dot(V1 * V2 * V3, E3))
@test isapprox(dot(v1 * v2 * v3, f1), dot(V1 * V2 * V3, F1))
@test isapprox(embed(v1 * v2 * v3), V1 * V2 * V3)
@test isapprox(embed(v1 * v2 * v3 * v4), V1 * V2 * V3 * V4)
@test isapprox(embed(exp(v1 * v2)), exp(V1 * V2))
@test isapprox(embed(bivector_exp(v1 * v2)), bivector_exp(V1 * V2))

run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)