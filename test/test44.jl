#Stand-alone tests of GA(4,4)
#All other algebras are compared to this.

v1 = inject(rand(8), GA.bas44)
v2 = inject(rand(8), GA.bas44)
v3 = inject(rand(8), GA.bas44)
v4 = inject(rand(8), GA.bas44)
v5 = inject(rand(8), GA.bas44)
v6 = inject(rand(8), GA.bas44)
v7 = inject(rand(8), GA.bas44)
v8 = inject(rand(8), GA.bas44)
v9 = inject(rand(8), GA.bas44)
v10 = inject(rand(8), GA.bas44)

me1 = v1*v2*v3*v4*v5*v6*v7*v8
me2 = v1*v2*v3*v4*v5*v6*v9*v10
me3 = v10*v9*v3*v4*v5*v6*v7*v8
mo1 = v1*v2*v3*v4*v5*v6*v7
mo2 = v1*v2*v3*v4*v8*v9*v10
mo3 = v8*v9*v10*v4*v5*v6*v7

#Distributivity
@test isapprox(me1*(me2+me3), me1*me2 + me1*me3)
@test isapprox(mo1*(me2+me3), mo1*me2 + mo1*me3)
@test isapprox(me1*(mo2+mo3), me1*mo2 + me1*mo3)
@test isapprox(mo1*(mo2+mo3), mo1*mo2 + mo1*mo3)

#Associativity
@test isapprox(me1*(me2*me3) , (me1*me2)*me3)
@test isapprox(mo1*(me2*me3) , (mo1*me2)*me3)
@test isapprox(me1*(mo2*me3) , (me1*mo2)*me3)
@test isapprox(me1*(me2*mo3) , (me1*me2)*mo3)
@test isapprox(mo1*(mo2*me3) , (mo1*mo2)*me3)
@test isapprox(mo1*(me2*mo3) , (mo1*me2)*mo3)
@test isapprox(me1*(mo2*mo3) , (me1*mo2)*mo3)
@test isapprox(mo1*(mo2*mo3) , (mo1*mo2)*mo3)


#Projection
@test isapprox(me1  , project(me1,0) + project(me1,2) + project(me1,4) + project(me1,6) + project(me1,8))
@test isapprox(mo1  , project(mo1,1) + project(mo1,3) + project(mo1,5) + project(mo1,7))


#Rotation
R = expb(v1*v2)
no1 = R*mo1*R'
no2 = R*mo2*R'
@test isapprox(dot(mo1,mo2),dot(no1,no2))