#Test suite for GA(3,3).
#Test stand-alone results and compares with GA(4,4)


e1 = GA.bas33[1]
e2 = GA.bas33[2]
e3 = GA.bas33[3]
f1 = GA.bas33[4]
f2 = GA.bas33[5]
f3 = GA.bas33[6]

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


#Comparison with GA(4,4)
bas44 = basis("GA44")
E1 = bas44[1]
E2 = bas44[2]
E3 = bas44[3]
F1 = bas44[5]
F2 = bas44[6]
F3 = bas44[7]
arr1 = rand(6)
v1 = inject(arr1,GA.bas33)
V1 = inject(arr1,[E1,E2,E3,F1,F2,F3])
arr2 = rand(6)
v2 = inject(arr2,GA.bas33)
V2 = inject(arr2,[E1,E2,E3,F1,F2,F3])
arr3 = rand(6)
v3 = inject(arr3,GA.bas33)
V3 = inject(arr3,[E1,E2,E3,F1,F2,F3])
arr4 = rand(6)
v4 = inject(arr4,GA.bas33)
V4 = inject(arr4,[E1,E2,E3,F1,F2,F3])
@test isapprox(dot(v1,v1),dot(V1,V1))
@test isapprox(dot(v1*v2*v3,e1), dot(V1*V2*V3,E1))
@test isapprox(dot(v1*v2*v3,e2), dot(V1*V2*V3,E2))
@test isapprox(dot(v1*v2*v3,e3), dot(V1*V2*V3,E3))
@test isapprox(dot(v1*v2*v3,f1), dot(V1*V2*V3,F1))
@test isapprox(embed(v1*v2*v3), V1*V2*V3)
@test isapprox(embed(v1*v2*v3*v4), V1*V2*V3*V4)
@test isapprox(embed(exp(v1*v2)),exp(V1*V2))
@test isapprox(embed(expb(v1*v2)),expb(V1*V2))


include("testcommon.jl")