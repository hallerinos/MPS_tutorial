# quickly compute nearest neighbor correlations & average magnetization
function measure_obs(AC)
    # define spin operators
    S0 = 1.0 .* [1    0;   0  1]
    Sx = 0.5 .* [0    1;   1  0]
    Sy = 0.5 .* [0 -1im; 1im  0]
    Sz = 0.5 .* [1    0;   0 -1]

    # nearest neighbor correlations
    @tensor xx[] := Sx[]*Sx[]
    @tensor yy[] := Sy[]*Sy[]
    @tensor zz[] := Sz[]*Sz[]
    xx_exp = @tensor AC[]*xx[]*conj(AC[])
    yy_exp = @tensor AC[]*yy[]*conj(AC[])
    zz_exp = @tensor AC[]*zz[]*conj(AC[])

    # average magnetization
    @tensor mx[] := 0.5*Sx[]*S0[] + 0.5*S0[]*Sx[]
    @tensor my[] := 0.5*Sy[]*S0[] + 0.5*S0[]*Sy[]
    @tensor mz[] := 0.5*Sz[]*S0[] + 0.5*S0[]*Sz[]
    mx_exp = @tensor AC[]*mx[]*conj(AC[])
    my_exp = @tensor AC[]*my[]*conj(AC[])
    mz_exp = @tensor AC[]*mz[]*conj(AC[])

    # return the two observables in two arrays
    return [xx_exp, yy_exp, zz_exp], [mx_exp, my_exp, mz_exp]
end