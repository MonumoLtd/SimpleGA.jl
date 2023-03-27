#Test suite for PGA.
#Test stand-alone results and compares with GA(4,4)

basPGA = PGA.basis
(e1, e2, e3, e0) = (basPGA[1], basPGA[2], basPGA[3], basPGA[4])
@test map(x -> dot(x, x), basPGA) == [1, 1, 1, 0]
@test testbas(basPGA)

I3 = e1 * e2 * e3
#! format:off
me1 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
me2 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
me3 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
mo1 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3
mo2 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3
mo3 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3
#! format:on

#Comparison with GA(4,4)
bas44 = GA44.basis
E1 = bas44[1]
E2 = bas44[2]
E3 = bas44[3]
E0 = bas44[4] + bas44[8]
arr1 = rand(4)
v1 = inject(arr1, basPGA)
V1 = inject(arr1, [E1, E2, E3, E0])
arr2 = rand(4)
v2 = inject(arr2, basPGA)
V2 = inject(arr2, [E1, E2, E3, E0])
arr3 = rand(4)
v3 = inject(arr3, basPGA)
V3 = inject(arr3, [E1, E2, E3, E0])
arr4 = rand(4)
v4 = inject(arr4, basPGA)
V4 = inject(arr4, [E1, E2, E3, E0])
@test isapprox(dot(v1, v1), dot(V1, V1))
@test isapprox(dot(v1 * v2 * v3, e1), dot(V1 * V2 * V3, E1))
@test isapprox(dot(v1 * v2 * v3, e2), dot(V1 * V2 * V3, E2))
@test isapprox(dot(v1 * v2 * v3, e3), dot(V1 * V2 * V3, E3))
@test isapprox(embed(v1 * v2 * v3), V1 * V2 * V3)
@test isapprox(embed(v1 * v2 * v3 * v4), V1 * V2 * V3 * V4)
@test isapprox(embed(exp(v1 * v2)), exp(V1 * V2))
@test isapprox(embed(bivector_exp(v1 * v2)), bivector_exp(V1 * V2))

run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)
