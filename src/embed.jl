#Routines for converting from small algebras to G(4,4)

function embed(a::GA20.Odd)
    bas44 = GA44.basis
    return real(a.c1) * bas44[1] + imag(a.c1) * bas44[2]
end

function embed(a::GA20.Even)
    bas44 = GA44.basis
    return real(a.c1) + imag(a.c1) * bas44[1] * bas44[2]
end

function embed(a::GA30.Even)
    bas44 = GA44.basis
    return a.w - a.x * bas44[2] * bas44[3] - a.y * bas44[3] * bas44[1] -
           a.z * bas44[1] * bas44[2]
end

function embed(a::GA30.Odd)
    return a.x * bas44[1] + a.y * bas44[2] + a.z * bas44[3] + a.w * bas[1] * bas[2] * bas[3]
end

function embed(a::STA.Even)
    bassta = STA.basis
    (g0, g1, g2, g3) = (bassta[1], bassta[2], bassta[3], bassta[4])
    bas44 = GA44.basis
    (G0, G1, G2, G3) = (bas44[1], bas44[5], bas44[6], bas44[7])
    i4 = g0 * g1 * g2 * g3
    I4 = G0 * G1 * G2 * G3
    evenlist = [
        tr(a),
        dot(a, g0 * g1),
        dot(a, g0 * g2),
        dot(a, g0 * g3),
        -dot(a, i4 * g0 * g1),
        -dot(a, i4 * g0 * g2),
        -dot(a, i4 * g0 * g3),
        -dot(a, i4),
    ]
    bas = [1, G0 * G1, G0 * G2, G0 * G3, I4 * G0 * G1, I4 * G0 * G2, I4 * G0 * G3, I4]
    return inject(evenlist, bas)
end

function embed(a::STA.Odd)
    bassta = STA.basis
    (g0, g1, g2, g3) = (bassta[1], bassta[2], bassta[3], bassta[4])
    bas44 = GA44.basis
    (G0, G1, G2, G3) = (bas44[1], bas44[5], bas44[6], bas44[7])
    i4 = g0 * g1 * g2 * g3
    I4 = G0 * G1 * G2 * G3
    oddlist = [
        dot(a, g0),
        -dot(a, g1),
        -dot(a, g2),
        -dot(a, g3),
        dot(a, i4 * g0),
        -dot(a, i4 * g1),
        -dot(a, i4 * g2),
        -dot(a, i4 * g3),
    ]
    bas = [G0, G1, G2, G3, I4 * G0, I4 * G1, I4 * G2, I4 * G3]
    return inject(oddlist, bas)
end

function embed(a::GA40.Even)
    bas40 = GA40.basis
    (e1, e2, e3, e4) = (bas40[1], bas40[2], bas40[3], bas40[4])
    bas44 = GA44.basis
    (E1, E2, E3, E4) = (bas44[1], bas44[2], bas44[3], bas44[4])
    b40 = [-e1 * e2, -e1 * e3, -e1 * e4, -e2 * e3, -e2 * e4, -e3 * e4, e1 * e2 * e3 * e4]
    b44 = [E1 * E2, E1 * E3, E1 * E4, E2 * E3, E2 * E4, E3 * E4, E1 * E2 * E3 * E4]
    res = tr(a)
    for i in 1:7
        res += dot(a, b40[i]) * b44[i]
    end
    return res
end

function embed(a::GA40.Odd)
    bas40 = GA40.basis
    (e1, e2, e3, e4) = (bas40[1], bas40[2], bas40[3], bas40[4])
    bas44 = GA44.basis
    (E1, E2, E3, E4) = (bas44[1], bas44[2], bas44[3], bas44[4])
    b40 = [e1, e2, e3, e4, -e2 * e3 * e4, -e1 * e3 * e4, -e1 * e2 * e4, -e1 * e2 * e3]
    b44 = [E1, E2, E3, E4, E2 * E3 * E4, E1 * E3 * E4, E1 * E2 * E4, E1 * E2 * E3]
    res = 0.0 * E1
    for i in 1:8
        res += dot(a, b40[i]) * b44[i]
    end
    return res
end

function embed(a::PGA.Even)
    bas44 = GA44.basis
    (E1, E2, E3, E0) = (bas44[1], bas44[2], bas44[3], bas44[4] + bas44[8])
    I3 = E1 * E2 * E3
    evenlist = [a.q.w, -a.q.z, -a.q.x, -a.q.y, -a.n.x, -a.n.y, -a.n.z, -a.n.w]
    bas = [1, E1 * E2, E2 * E3, E3 * E1, E0 * E1, E0 * E2, E0 * E3, E0 * I3]
    return inject(evenlist, bas)
end

function embed(a::PGA.Odd)
    bas44 = GA44.basis
    (E1, E2, E3, E0) = (bas44[1], bas44[2], bas44[3], bas44[4] + bas44[8])
    I3 = E1 * E2 * E3
    oddlist = [a.n.w, -a.q.x, -a.q.y, -a.q.z, a.n.z, a.n.y, a.n.x, -a.q.w]
    bas = [E0, E1, E2, E3, E0 * E2 * E1, E0 * E1 * E3, E0 * E3 * E2, I3]
    return inject(oddlist, bas)
end

function embed(a::CGA.Even)
    bascga = CGA.basis
    (e1, e2, e3, e4, f4) = (bascga[1], bascga[2], bascga[3], bascga[4], bascga[5])
    bas44 = GA44.basis
    (E1, E2, E3, E4, F4) = (bas44[1], bas44[2], bas44[3], bas44[4], bas44[8])
    i5 = e1 * e2 * e3 * e4 * f4
    I5 = E1 * E2 * E3 * E4 * F4
    #! format:off
    evenbas = [-e1*e2, -e1*e3, -e1*e4, -e2*e3, -e2*e4, -e3*e4, e1*f4, e2*f4, e3*f4, e4*f4,
                -i5*e1, -i5*e2, -i5*e3, -i5*e4, i5*f4]
    #! format:on
    vals = map(bs -> dot(a, bs), evenbas)
    #! format:off
    resbas = [E1*E2, E1*E3, E1*E4, E2*E3, E2*E4, E3*E4, E1*F4, E2*F4, E3*F4, E4*F4,
                I5*E1, I5*E2, I5*E3, I5*E4, I5*F4]
    #! format:on
    return tr(a) + inject(vals, resbas)
end

function embed(a::CGA.Odd)
    bascga = CGA.basis
    (e1, e2, e3, e4, f4) = (bascga[1], bascga[2], bascga[3], bascga[4], bascga[5])
    bas44 = GA44.basis
    (E1, E2, E3, E4, F4) = (bas44[1], bas44[2], bas44[3], bas44[4], bas44[8])
    i5 = e1 * e2 * e3 * e4 * f4
    I5 = E1 * E2 * E3 * E4 * F4
    #! format:off
    oddbas = [e1, e2, e3, e4, -f4,
            i5*e1*e2, i5*e1*e3, i5*e1*e4, i5*e2*e3, i5*e2*e4, i5*e3*e4, -i5*e1*f4, -i5*e2*f4, -i5*e3*f4, -i5*e4*f4,
            -i5 ]
    #! format:on
    vals = map(bs -> dot(a, bs), oddbas)
    #! format:off
    resbas = [E1, E2, E3, E4, F4,
                I5*E1*E2, I5*E1*E3, I5*E1*E4, I5*E2*E3, I5*E2*E4, I5*E3*E4, I5*E1*F4, I5*E2*F4, I5*E3*F4, I5*E4*F4,
                I5]
    #! format:on
    return inject(vals, resbas)
end

function embed(a::GA33.Even)
    bas33 = GA33.basis
    (e1, e2, e3, f1, f2, f3) = (bas33[1], bas33[2], bas33[3], bas33[4], bas33[5], bas33[6])
    i6 = e1 * e2 * e3 * f1 * f2 * f3
    #! format:off
    evenbas = [ -e1*e2, -e1*e3, -e2*e3, -f1*f2, -f1*f3, -f2*f3,
                e1*f1, e1*f2, e1*f3, e2*f1, e2*f2, e2*f3, e3*f1, e3*f2, e3*f3,
                -e1*e2*i6, -e1*e3*i6, -e2*e3*i6, -f1*f2*i6, -f1*f3*i6, -f2*f3*i6,
                e1*f1*i6, e1*f2*i6, e1*f3*i6,e2*f1*i6, e2*f2*i6, e2*f3*i6, e3*f1*i6, e3*f2*i6, e3*f3*i6,
                i6]
        #! format:on
    vals = map(bs -> dot(a, bs), evenbas)
    bas44 = GA44.basis
    (E1, E2, E3, F1, F2, F3) = (bas44[1], bas44[2], bas44[3], bas44[5], bas44[6], bas44[7])
    I6 = E1 * E2 * E3 * F1 * F2 * F3
    #! format:off
    resbas = [  E1*E2, E1*E3, E2*E3, F1*F2, F1*F3, F2*F3,
                E1*F1, E1*F2, E1*F3, E2*F1, E2*F2, E2*F3, E3*F1, E3*F2, E3*F3,
                E1*E2*I6, E1*E3*I6, E2*E3*I6, F1*F2*I6, F1*F3*I6, F2*F3*I6,
                E1*F1*I6, E1*F2*I6, E1*F3*I6,E2*F1*I6, E2*F2*I6, E2*F3*I6, E3*F1*I6, E3*F2*I6, E3*F3*I6,
                I6]
    #! format:on
    return tr(a) + inject(vals, resbas)
end

function embed(a::GA33.Odd)
    bas33 = GA33.basis
    (e1, e2, e3, f1, f2, f3) = (bas33[1], bas33[2], bas33[3], bas33[4], bas33[5], bas33[6])
    i6 = e1 * e2 * e3 * f1 * f2 * f3
    #! format:off
    oddbas = [  e1, e2, e3, -f1, -f2, -f3,
                -e1*e2*e3,
                e1*e2*f1, e1*e3*f1, e2*e3*f1, e1*e2*f2, e1*e3*f2, e2*e3*f2, e1*e2*f3, e1*e3*f3, e2*e3*f3,
                -e1*f1*f2, - e1*f1*f3, -e1*f2*f3, -e2*f1*f2, - e2*f1*f3, -e2*f2*f3, -e3*f1*f2, - e3*f1*f3, -e3*f2*f3,
                f1*f2*f3,
                -e1*i6, -e2*i6, -e3*i6, f1*i6, f2*i6, f3*i6]
    #! format:on
    vals = map(bs -> dot(a, bs), oddbas)
    bas44 = GA44.basis
    (E1, E2, E3, F1, F2, F3) = (bas44[1], bas44[2], bas44[3], bas44[5], bas44[6], bas44[7])
    I6 = E1 * E2 * E3 * F1 * F2 * F3
    #! format:off
    resbas = [  E1, E2, E3, F1, F2, F3,
                E1*E2*E3,
                E1*E2*F1, E1*E3*F1, E2*E3*F1, E1*E2*F2, E1*E3*F2, E2*E3*F2, E1*E2*F3, E1*E3*F3, E2*E3*F3,
                E1*F1*F2,  E1*F1*F3, E1*F2*F3, E2*F1*F2, E2*F1*F3, E2*F2*F3, E3*F1*F2,  E3*F1*F3, E3*F2*F3,
                F1*F2*F3,
                E1*I6, E2*I6, E3*I6, F1*I6, F2*I6, F3*I6]
    #! format:on
    return inject(vals, resbas)
end

function embed(a::GA24.Even)
    bas24 = GA24.basis
    (g0, g1, g2, g3, g4, g5) = (bas24[1], bas24[2], bas24[3], bas24[4], bas24[5], bas24[6])
    i6 = g0 * g1 * g2 * g3 * g4 * g5
    #! format:off
    evenbas = [ g1*g0, g2*g0, g3*g0, g4*g0, -g5*g0,
                -g1*g2, -g1*g3, -g1*g4, g1*g5,
                -g2*g3, -g2*g4, g2*g5,
                -g3*g4, g3*g5,
                g4*g5,
                -g1*g0*i6, -g2*g0*i6, -g3*g0*i6, -g4*g0*i6, g5*g0*i6,
                g1*g2*i6, g1*g3*i6, g1*g4*i6, -g1*g5*i6,
                g2*g3*i6, g2*g4*i6, -g2*g5*i6,
                g3*g4*i6, -g3*g5*i6,
                -g4*g5*i6,
                -i6]
    #! format:on
    vals = map(bs -> dot(a, bs), evenbas)
    bas44 = GA44.basis
    (G0, G1, G2, G3, G4, G5) = (bas44[1], bas44[5], bas44[6], bas44[7], bas44[8], bas44[4])
    I6 = G0 * G1 * G2 * G3 * G4 * G5
    #! format:off
    resbas = [  G1*G0, G2*G0, G3*G0, G4*G0, G5*G0,
                G1*G2, G1*G3, G1*G4, G1*G5,
                G2*G3, G2*G4, G2*G5,
                G3*G4, G3*G5,
                G4*G5,
                G1*G0*I6, G2*G0*I6, G3*G0*I6, G4*G0*I6, G5*G0*I6,
                G1*G2*I6, G1*G3*I6, G1*G4*I6, G1*G5*I6,
                G2*G3*I6, G2*G4*I6, G2*G5*I6,
                G3*G4*I6, G3*G5*I6,
                G4*G5*I6,
                I6 ]
   #! format:on
    return tr(a) + inject(vals, resbas)
end

function embed(a::GA24.Odd)
    bas24 = GA24.basis
    (g0, g1, g2, g3, g4, g5) = (bas24[1], bas24[2], bas24[3], bas24[4], bas24[5], bas24[6])
    i6 = g0 * g1 * g2 * g3 * g4 * g5
    #! format:off
    oddbas = [  g0, -g1, -g2, -g3, -g4, g5,
                -g0*g1*g2, -g0*g1*g3, -g0*g1*g4, g0*g1*g5, -g0*g2*g3, -g0*g2*g4, g0*g2*g5, -g0*g3*g4, g0*g3*g5, g0*g4*g5,
                -g0*g1*g2*i6, -g0*g1*g3*i6, -g0*g1*g4*i6, g0*g1*g5*i6, -g0*g2*g3*i6, -g0*g2*g4*i6, g0*g2*g5*i6, -g0*g3*g4*i6, g0*g3*g5*i6, g0*g4*g5*i6,
                g0*i6, -g1*i6, -g2*i6, -g3*i6, -g4*i6, g5*i6 ]
    #! format:off
    vals = map(bs->dot(a,bs),oddbas)
    bas44 = GA44.basis
    (G0,G1,G2,G3,G4,G5) = (bas44[1], bas44[5], bas44[6], bas44[7], bas44[8], bas44[4])
    I6 = G0*G1*G2*G3*G4*G5
    #! format:off
    resbas = [  G0, G1, G2, G3, G4, G5,
                G0*G1*G2, G0*G1*G3, G0*G1*G4, G0*G1*G5, G0*G2*G3, G0*G2*G4, G0*G2*G5, G0*G3*G4, G0*G3*G5, G0*G4*G5,
                G0*G1*G2*I6, G0*G1*G3*I6, G0*G1*G4*I6, G0*G1*G5*I6, G0*G2*G3*I6, G0*G2*G4*I6, G0*G2*G5*I6, G0*G3*G4*I6, G0*G3*G5*I6, G0*G4*G5*I6,
                G0*I6, G1*I6, G2*I6, G3*I6, G4*I6, G5*I6 ]
    #! format:on
    return inject(vals, resbas)
end
