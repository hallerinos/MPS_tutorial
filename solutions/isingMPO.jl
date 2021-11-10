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
    W = zeros(3,2,3,2)
    σ₁ = [+0.0 +1.0; +1.0 +0.0]
    σ₃ = [+1.0 +0.0; +0.0 -1.0]
    σ₀ = Matrix(I, 2, 2)

    W[1,:,1,:] = σ₀
    W[1,:,2,:] = σ₁
    W[1,:,3,:] = h*σ₃
    W[2,:,3,:] = J*σ₁
    W[3,:,3,:] = σ₀

    L = zeros(1,3)
    L[1] = 1.0
    R = zeros(3,1)
    R[end] = 1.0

    return L, W, R
end

function MPO_check()
    σ₁ = [+0.0 +1.0; +1.0 +0.0]
    σ₃ = [+1.0 +0.0; +0.0 -1.0]
    σ₀ = Matrix(I, 2, 2)

    J=1; h=2
    L, W, R = isingMPO(J=J, h=h)

    ⊗(A,B) = kron(A,B)
    @tensor LWWR[:] := L[4,1]*W[1,-1,2,-3]*W[2,-2,3,-4]*R[3,4]
    @show J*(σ₁⊗σ₁)+h*(σ₃⊗σ₀+σ₀⊗σ₃) ≈ reshape(LWWR, (4,4))

    @tensor LWWWR[:] := L[5,1]*W[1,-4,2,-1]*W[2,-5,3,-2]*W[3,-6,4,-3]*R[4,5]
    @show J*(σ₁⊗σ₁⊗σ₀+σ₀⊗σ₁⊗σ₁)+h*(σ₃⊗σ₀⊗σ₀+σ₀⊗σ₃⊗σ₀+σ₀⊗σ₀⊗σ₃) ≈ reshape(LWWWR, (8,8))

    @tensor LWWWWR[:] := L[6,1]*W[1,-5,2,-1]*W[2,-6,3,-2]*W[3,-7,4,-3]*W[4,-8,5,-4]*R[5,6]
    @show J*(σ₁⊗σ₁⊗σ₀⊗σ₀+σ₀⊗σ₁⊗σ₁⊗σ₀+σ₀⊗σ₀⊗σ₁⊗σ₁)+h*(σ₃⊗σ₀⊗σ₀⊗σ₀+σ₀⊗σ₃⊗σ₀⊗σ₀+σ₀⊗σ₀⊗σ₃⊗σ₀+σ₀⊗σ₀⊗σ₀⊗σ₃) ≈ reshape(LWWWWR, (16,16))
    return 
end

# MPO_check()