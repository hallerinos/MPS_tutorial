include("isingMPO.jl")
include("contractions.jl")
using KrylovKit  # we use "eigsolve" from this package

function random_guess(χ, d)
    AC = randn(χ, d, d, χ); AC /= norm(AC)
end

function smart_guess(U, Λ, Λ_prev, Vdag)
    @tensor AC[:] := Λ[-1,1]*Vdag[1,-2,2]*pinv(Λ_prev)[2,3]*U[3,-3,4]*Λ[4,-4]
end

function iDMRG(;verbose=1, J=1, h=1, n=100, χ=16, kwargs...)
    # get the MPO
    L, W, R = isingMPO(J=J, h=h)
    LW = reshape(L, (1,3,1))
    RW = reshape(R, (1,3,1))

    χis = 1
    Λ = ones(1,1)

    ene = 0.0
    U = randn(1, size(W,4), 1)
    Vdag = randn(1, size(W,4), 1)
    Λ_prev = Λ
    AC = randn(1, size(W,4), size(W,4), 1)
    Δρ = 0.0
    @time for iteration=1:n
        # check energy change wrt previous iteration
        ene_prev = ene

        # prepare guess
        # AC = random_guess(χis, size(W,4))
        AC = smart_guess(U, Λ, Λ_prev, Vdag)

        # solve for the lowest eigenvalue & eigenvector
        eval, evec = eigsolve(AC, 1, :SR) do x grad(x, W, LW, RW) end

        # extract it from the list of return values
        ene = eval[1]
        AC = evec[1]

        # energy density
        ene_per_site = (ene-ene_prev)/2

        # save previous Λ matrix
        Λ_prev = Λ

        # perform SVD
        ACmat = reshape(AC, (prod(size(AC)[1:2]), prod(size(AC)[3:end])))  # reshape to matrix
        U, Λ, V = svd(ACmat)  # extact the left and right isometries
        Sfull = -sum(Λ.^2 .* log.(Λ.^2))  # von Neumann entanglement entropy before truncation
        χis = min(χ,length(Λ))  # check if number of singular values is less than χ
        U = U[:,1:χis]  # truncate left isometry
        Δρ = sum(Λ[χis+1:end].^2)  # "truncated weight", it gives an estimate for the overall error
        Λ = Λ[1:χis]  # truncate singular values
        Λ /= norm(Λ)  # restore norm
        Strunc = -sum(Λ.^2 .* log.(Λ.^2))  # von Neumann entanglement entropy after truncation
        Λ = Diagonal(Λ)  # cast it to a diagonal matrix
        Vdag = V'[1:χis,:]  # truncate right isometry

        ΔAC = sqrt(norm(ACmat .- U*Λ*Vdag))  # optimal state vs compressed state
        ΔS = Sfull - Strunc  # the entropy is well approximated

        if verbose==1  # to make output less annoying later on
            @show iteration ene_per_site  # display iteration step and energy density
            @show χis Δρ ΔAC ΔS  # show current bond dimension and some error measures
        end

        # need to reshape U and Vdag into tensor form to continue
        U = reshape(U, (size(AC,1), size(AC,2), χis))
        Vdag = reshape(Vdag, (χis, size(AC,3), size(AC,4)))

        # update environments
        LW = absorbU(U, W, LW)
        RW = absorbVdag(Vdag, W, RW)
    end

    return AC, real(ene)
end

# iDMRG()