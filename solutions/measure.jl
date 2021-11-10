# quickly compute nearest neighbor correlations & average magnetization
function measure_obs(AC)
    # define spin operators
    S0 = 1.0 .* [1    0;   0  1]
    Sx = 0.5 .* [0    1;   1  0]
    Sy = 0.5 .* [0 -1im; 1im  0]
    Sz = 0.5 .* [1    0;   0 -1]

    # nearest neighbor correlations
    @tensor xx[-1,-2,-3,-4] := Sx[-1,-3]*Sx[-2,-4]
    xx_exp = @tensor AC[1 2 3 4]*xx[2 3 5 6]*conj(AC[1 5 6 4])

    # average magnetization
    @tensor mz[-1,-2,-3,-4] := 0.5*Sz[-1,-3]*S0[-2,-4] + 0.5*S0[-1,-3]*Sz[-2,-4]
    mz_exp = @tensor AC[1 2 3 4]*mz[2 3 5 6]*conj(AC[1 5 6 4])

    # return the two observables in two arrays
    return xx_exp, mz_exp
end