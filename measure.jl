# quickly compute nearest neighbor correlations & average magnetization
function measure_obs(AC)
    # define spin operators
    S0 = 1.0 .* [1    0;   0  1]
    Sx = 0.5 .* [0    1;   1  0]
    Sy = 0.5 .* [0 -1im; 1im  0]
    Sz = 0.5 .* [1    0;   0 -1]

    # nearest neighbor correlations
    @tensor xx[-1,-2,-3,-4] := Sx[-1,-3]*Sx[-2,-4]
    @tensor yy[-1,-2,-3,-4] := Sy[-1,-3]*Sy[-2,-4]
    @tensor zz[-1,-2,-3,-4] := Sz[-1,-3]*Sz[-2,-4]
    xx_exp = @tensor AC[1 2 3 4]*xx[2 3 5 6]*conj(AC[1 5 6 4])
    yy_exp = @tensor AC[1 2 3 4]*yy[2 3 5 6]*conj(AC[1 5 6 4])
    zz_exp = @tensor AC[1 2 3 4]*zz[2 3 5 6]*conj(AC[1 5 6 4])

    # average magnetization
    @tensor mx[-1,-2,-3,-4] := 0.5*Sx[-1,-3]*S0[-2,-4] + 0.5*S0[-1,-3]*Sx[-2,-4]
    @tensor my[-1,-2,-3,-4] := 0.5*Sy[-1,-3]*S0[-2,-4] + 0.5*S0[-1,-3]*Sy[-2,-4]
    @tensor mz[-1,-2,-3,-4] := 0.5*Sz[-1,-3]*S0[-2,-4] + 0.5*S0[-1,-3]*Sz[-2,-4]
    mx_exp = @tensor AC[1 2 3 4]*mx[2 3 5 6]*conj(AC[1 5 6 4])
    my_exp = @tensor AC[1 2 3 4]*my[2 3 5 6]*conj(AC[1 5 6 4])
    mz_exp = @tensor AC[1 2 3 4]*mz[2 3 5 6]*conj(AC[1 5 6 4])

    # return the two observables in two arrays
    return [xx_exp, yy_exp, zz_exp], [mx_exp, my_exp, mz_exp]
end