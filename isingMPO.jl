using TensorOperations, LinearAlgebra

function isingMPO(;J=1.0, h=1.0, kwargs...)
    # index order is fixed to
    #       -4       
    #        |       
    #       ---      
    # -1 --| W |-- -3
    #       ---      
    #        |       
    #       -2       
    σ₁ = [+0.0 +1.0; +1.0 +0.0]
    σ₃ = [+1.0 +0.0; +0.0 -1.0]
    σ₀ = Matrix(I, 2, 2)

    # W = ?
    # L = ?
    # R = ?

    return L, W, R
end

function MPO_check()
    σ₁ = [+0.0 +1.0; +1.0 +0.0]
    σ₃ = [+1.0 +0.0; +0.0 -1.0]
    σ₀ = Matrix(I, 2, 2)

    J=1; h=2
    # L, W, R = isingMPO(J=J, h=h)  # de-comment this line after you defined W, L, R

    ⊗(A,B) = kron(A,B)
    # # now for the sanity checks -- implement the contractions and check by de-commenting line after line
    # @tensor LWWR[:] := L[]*W[]*W[]*R[]
    # @show J*(σ₁⊗σ₁)+h*(σ₃⊗σ₀+σ₀⊗σ₃) ≈ reshape(LWWR, (4,4))

    # @tensor LWWWR[:] := L[]*W[]*W[]*W[]*R[]
    # @show J*(σ₁⊗σ₁⊗σ₀+σ₀⊗σ₁⊗σ₁)+h*(σ₃⊗σ₀⊗σ₀+σ₀⊗σ₃⊗σ₀+σ₀⊗σ₀⊗σ₃) ≈ reshape(LWWWR, (8,8))

    # @tensor LWWWWR[:] := L[]*W[]*W[]*W[]*W[]*R[]
    # @show J*(σ₁⊗σ₁⊗σ₀⊗σ₀+σ₀⊗σ₁⊗σ₁⊗σ₀+σ₀⊗σ₀⊗σ₁⊗σ₁)+h*(σ₃⊗σ₀⊗σ₀⊗σ₀+σ₀⊗σ₃⊗σ₀⊗σ₀+σ₀⊗σ₀⊗σ₃⊗σ₀+σ₀⊗σ₀⊗σ₀⊗σ₃) ≈ reshape(LWWWWR, (16,16))
    return 
end

# sanity check if construction of MPO is correct (comment next line if you finished the task)
MPO_check()