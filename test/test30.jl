#Test suite for GA(3,0).
#Test stand-alone results and compares with GA(4,4)

bas30 = GA30.basis
e1 = bas30[1]
e2 = bas30[2]
e3 = bas30[3]

@test map(x -> dot(x, x), bas30) == [1, 1, 1]
@test testbas(bas30)

me1 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
me2 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
me3 = rand() + rand() * e1 * e2 + e1 * e3 * rand() + e3 * rand() * e2
mo1 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1
mo2 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1
mo3 = rand() * e1 + rand() * e2 + e3 * rand() + e3 * rand() * e2 * e1

#Comparison with GA(4,4)
bas44 = GA44.basis
(E1, E2, E3) = (bas44[1], bas44[2], bas44[3])
arr1 = rand(3)
v1 = inject(arr1, bas30)
V1 = inject(arr1, [E1, E2, E3])
arr2 = rand(3)
v2 = inject(arr2, bas30)
V2 = inject(arr2, [E1, E2, E3])
arr3 = rand(3)
v3 = inject(arr3, bas30)
V3 = inject(arr3, [E1, E2, E3])
@test isapprox(dot(v1, v1), dot(V1, V1))
@test isapprox(dot(v1 * v2 * v3, e1), dot(V1 * V2 * V3, E1))
@test isapprox(dot(v1 * v2 * v3, e2), dot(V1 * V2 * V3, E2))
@test isapprox(dot(v1 * v2 * v3, e3), dot(V1 * V2 * V3, E3))
@test isapprox(embed(exp(v1 * v2)), exp(V1 * V2))
@test isapprox(embed(bivector_exp(v1 * v2)), bivector_exp(V1 * V2))

run_common_tests(me1, me2, me3, mo1, mo2, mo3, v1, v2)

# Conversion
run_conversion_tests(me1, me2, mo1, mo2, GA30.Even{Float32}, GA30.Odd{Float32})