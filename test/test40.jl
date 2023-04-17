#Test suite for GA(4,0).
#Test stand-alone results and compares with GA(4,4)

bas40 = GA40.basis
(e1, e2, e3, e4) = (bas40[1], bas40[2], bas40[3], bas40[4])
E4 = e1 * e2 * e3 * e4

@test map(x -> dot(x, x), bas40) == [1, 1, 1, 1]
@test testbas(bas40)

#! format:off
me1 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + rand()*e3*e4 + rand()*E4
me2 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + rand()*e3*e4 + rand()*E4
me3 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + rand()*e3*e4 + E4/rand()
mo1 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand() + E4*(rand()*e1 + rand()*e2 + e3*rand() + e4*rand())
mo2 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand() + E4*(rand()*e1 + rand()*e2 + e3*rand() + e4*rand())
mo3 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand() + E4*(rand()*e1 + rand()*e2 + e3*rand() + e4*rand())
#! format:on

#Comparison with GA(4,4)
bas44 = GA44.basis
(E1, E2, E3, E4) = (bas44[1], bas44[2], bas44[3], bas44[4])
arr1 = rand(4)
v1 = inject(arr1, bas40)
V1 = inject(arr1, [E1, E2, E3, E4])
arr2 = rand(4)
v2 = inject(arr2, bas40)
V2 = inject(arr2, [E1, E2, E3, E4])
arr3 = rand(4)
v3 = inject(arr3, bas40)
V3 = inject(arr3, [E1, E2, E3, E4])
arr4 = rand(4)
v4 = inject(arr4, bas40)
V4 = inject(arr4, [E1, E2, E3, E4])
@test isapprox(dot(v1, v1), dot(V1, V1))
@test isapprox(dot(v1 * v2 * v3, e1), dot(V1 * V2 * V3, E1))
@test isapprox(dot(v1 * v2 * v3, e2), dot(V1 * V2 * V3, E2))
@test isapprox(dot(v1 * v2 * v3, e3), dot(V1 * V2 * V3, E3))
@test isapprox(embed(v1 * v2 * v3), V1 * V2 * V3)
@test isapprox(embed(v1 * v2 * v3 * v4), V1 * V2 * V3 * V4)
@test isapprox(embed(exp(v1 * v2)), exp(V1 * V2))
@test isapprox(embed(bivector_exp(v1 * v2)), bivector_exp(V1 * V2))

run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)

# Conversion
run_conversion_tests(me1, me2, mo1, mo2, GA40.Even{Float32}, GA40.Odd{Float32})
