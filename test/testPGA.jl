#Test suite for PGA.
#Test stand-alone results and compares with GA(4,4)

baspga = basis("PGA")
e1 = baspga[1]
e2 = baspga[2]
e3 = baspga[3]
e0 = baspga[4]

I3 = e1*e2*e3

me1 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
me2 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
me3 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e0 + rand()*e2*e0 + rand()*e3*e0 + rand()*I3*e0
mo1 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3
mo2 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3
mo3 = rand()*e1 + rand()*e2 + e3*rand() + e0/rand() + I3*e0*(rand()*e1 + rand()*e2 + e3*rand()) + + rand()*I3


#Comparison with GA(4,4)
bas44 = basis("GA44")
E1 = bas44[1]
E2 = bas44[2]
E3 = bas44[3]
E0 = bas44[4] + bas44[8]
arr1 = rand(4)
v1 = inject(arr1,baspga)
V1 = inject(arr1,[E1,E2,E3,E0])
arr2 = rand(4)
v2 = inject(arr2,baspga)
V2 = inject(arr2,[E1,E2,E3,E0])
arr3 = rand(4)
v3 = inject(arr3,baspga)
V3 = inject(arr3,[E1,E2,E3,E0])
arr4 = rand(4)
v4 = inject(arr4,baspga)
V4 = inject(arr4,[E1,E2,E3,E0])
@test isapprox(dot(v1,v1),dot(V1,V1))
@test isapprox(dot(v1*v2*v3,e1), dot(V1*V2*V3,E1))
@test isapprox(dot(v1*v2*v3,e2), dot(V1*V2*V3,E2))
@test isapprox(dot(v1*v2*v3,e3), dot(V1*V2*V3,E3))
@test isapprox(embed(v1*v2*v3), V1*V2*V3)
@test isapprox(embed(v1*v2*v3*v4), V1*V2*V3*V4)
@test isapprox(embed(exp(v1*v2)),exp(V1*V2))
@test isapprox(embed(expb(v1*v2)),expb(V1*V2))

include("testcommon.jl")