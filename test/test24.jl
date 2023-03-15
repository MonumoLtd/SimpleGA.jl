#Test suite for GA(2,4).
#Test stand-alone results and compares with GA(2,4)


g0 = GA.bas24[1]
g1 = GA.bas24[2]
g2 = GA.bas24[3]
g3 = GA.bas24[4]
g4 = GA.bas24[5]
g5 = GA.bas24[6]


me1 = rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4 +
    g2*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*g5*(rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4)
me2 = rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4 +
    g2*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*g5*(rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4)
me3 = rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4 +
    g2*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*g5*(rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4)
mo1 = (rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g5*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g2*g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ))
mo2 = (rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g5*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g2*g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ))
mo3 = (rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ) ) +
    g2*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g5*( rand() + rand()*g0*g1 + g0*g3*rand() + g0*rand()*g4 - rand()*g3*g1 + g1*g4/rand() + rand()*g3*g4 + rand()*g0*g1*g3*g4) +
    g2*g5*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() + g0*g1*g3*g4*(rand()*g0 + rand()*g1 + g3*rand() + g4*rand() ))




#Comparison with GA(4,4)
bas44 = basis("GA44")
G0 = bas44[1]
G1 = bas44[5]
G2 = bas44[6]
G3 = bas44[7]
G4 = bas44[8]
G5 = bas44[4]
arr1 = rand(6)
v1 = inject(arr1,GA.bas24)
V1 = inject(arr1,[G0,G1,G2,G3,G4,G5])
arr2 = rand(6)
v2 = inject(arr2,GA.bas24)
V2 = inject(arr2,[G0,G1,G2,G3,G4,G5])
arr3 = rand(6)
v3 = inject(arr3,GA.bas24)
V3 = inject(arr3,[G0,G1,G2,G3,G4,G5])
arr4 = rand(6)
v4 = inject(arr4,GA.bas24)
V4 = inject(arr4,[G0,G1,G2,G3,G4,G5])
@test isapprox(dot(v1,v1),dot(V1,V1))
@test isapprox(dot(v1*v2*v3,g1), dot(V1*V2*V3,G1))
@test isapprox(dot(v1*v2*v3,g2), dot(V1*V2*V3,G2))
@test isapprox(dot(v1*v2*v3,g3), dot(V1*V2*V3,G3))
@test isapprox(dot(v1*v2*v3,g0), dot(V1*V2*V3,G0))
@test isapprox(embed(v1*v2*v3), V1*V2*V3)
@test isapprox(embed(v1*v2*v3*v4), V1*V2*V3*V4)
@test isapprox(embed(exp(v1*v2)),exp(V1*V2))
@test isapprox(embed(expb(v1*v2)),expb(V1*V2))

    
include("testcommon.jl")
