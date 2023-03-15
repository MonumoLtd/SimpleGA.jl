#Test suite for CGA.
#Test stand-alone results and compares with GA(4,4)

e1= GA.basCGA[1]
e2= GA.basCGA[2]
e3= GA.basCGA[3]
e4= GA.basCGA[4]
f4= GA.basCGA[5]
I5 = e1*e2*e3*e4*f4

#CGA
me1 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + 
f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4) + I5*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4 + rand()*f4)
me2 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + 
f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4) + I5*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4 + rand()*f4)
me3 = rand() + rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + 
f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4) + I5*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4 + rand()*f4)
mo1 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand()+ f4/rand() + rand()*I5 +
I5*(rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4))
mo2 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand()+ f4/rand() + rand()*I5 +
I5*(rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4))
mo3 = rand()*e1 + rand()*e2 + e3*rand() + e4*rand()+ f4/rand() + rand()*I5 +
I5*(rand()*e1*e2 + e1*e3*rand() + e3*rand()*e2 - rand()*e1*e4 + rand()*e2*e4 + e3*e4/rand() + f4*(rand()*e1 + rand()*e2 +rand()*e3 +rand()*e4))


#Comparison with GA(4,4)
bas44 = basis("GA44")
E1 = bas44[1]
E2 = bas44[2]
E3 = bas44[3]
E4 = bas44[4]
F4 = bas44[8]
arr1 = rand(5)
v1 = inject(arr1,GA.basCGA)
V1 = inject(arr1,[E1,E2,E3,E4,F4])
arr2 = rand(5)
v2 = inject(arr2,GA.basCGA)
V2 = inject(arr2,[E1,E2,E3,E4,F4])
arr3 = rand(5)
v3 = inject(arr3,GA.basCGA)
V3 = inject(arr3,[E1,E2,E3,E4, F4])
arr4 = rand(5)
v4 = inject(arr4,GA.basCGA)
V4 = inject(arr4,[E1,E2,E3,E4,F4])
arr5 = rand(5)
v5 = inject(arr5,GA.basCGA)
V5 = inject(arr5,[E1,E2,E3,E4,F4])
@test isapprox(dot(v1,v1),dot(V1,V1))
@test isapprox(dot(v1*v2*v3,e1), dot(V1*V2*V3,E1))
@test isapprox(dot(v1*v2*v3,e2), dot(V1*V2*V3,E2))
@test isapprox(dot(v1*v2*v3,e3), dot(V1*V2*V3,E3))
@test isapprox(dot(v1*v2*v3,f4), dot(V1*V2*V3,F4))
@test isapprox(embed(v1*v2*v3), V1*V2*V3)
@test isapprox(embed(v1*v2*v3*v4), V1*V2*V3*V4)
@test isapprox(embed(v1*v2*v3*v4*v5), V1*V2*V3*V4*V5)
@test isapprox(embed(exp(v1*v2)),exp(V1*V2))
@test isapprox(embed(expb(v1*v2)),expb(V1*V2))

include("testcommon.jl")