#Test shared across all algebras.

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
@test isapprox(me1  , project(me1,0) + project(me1,2) + project(me1,4) + project(me1,6))
@test isapprox(mo1  , project(mo1,1) + project(mo1,3) + project(mo1,5))


#Rotation
R = expb(v1*v2)
no1 = R*mo1*R'
no2 = R*mo2*R'
@test isapprox(dot(mo1,mo2),dot(no1,no2))
ne1 = R*me1*R'
ne2 = R*me2*R'
@test isapprox(dot(me1,me2),dot(ne1,ne2))


#Reverse
@test isapprox((me1+me1')/2, tr(me1) + project(me1,4))
@test isapprox((me1-me1')/2, project(me1,2) + project(me1,6))
@test isapprox((mo1+mo1')/2, project(mo1,1) + project(mo1,5))
@test isapprox((mo1-mo1')/2, project(mo1,3))
