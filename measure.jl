# quickly compute nearest neighbor correlations & average magnetization
function measure_obs(AC)
    # define spin operators
    S0 = 1.0 .* [1    0;   0  1]
    Sx = 0.5 .* [0    1;   1  0]
    Sy = 0.5 .* [0 -1im; 1im  0]
    Sz = 0.5 .* [1    0;   0 -1]

    # nearest neighbor correlations
    @tensor xx[] := Sx[]*Sx[]
    xx_exp = @tensor AC[]*xx[]*conj(AC[])

    # average magnetization
    @tensor mz[] := 0.5*Sz[]*S0[] + 0.5*S0[]*Sz[]
    mz_exp = @tensor AC[]*mz[]*conj(AC[])

    # return the two observables in two arrays
    return xx_exp, mz_exp
end