using SparseArrays, LinearAlgebra
include("generate_spin_matrices.jl")

⊗(A,B) = kron(A, B)

function generate_spin_ops(s::Float64, N::Int64)
    spin_matrices = generate_spin_matrices(s)
    d = size(spin_matrices["Id"],1)

    Ô = Dict()

    ρ = spin_matrices["S+"]*spin_matrices["S-"]
    N̂ = spzeros(Float64, d^N, d^N)
    for n=1:N, sm in keys(spin_matrices)
        idL = sparse(I, d^(n-1), d^(n-1))
        idR = sparse(I, d^(N-n), d^(N-n))
        mat = idL ⊗ spin_matrices[sm] ⊗ idR
        Ô[sm, n] = mat

        N̂ += idL ⊗ ρ ⊗ idR
    end
    return Ô, N̂
end

function generate_string2(N::Int64)
    sz = -[1, -1]
    SZ = 1
    for n=1:N
        SZ = [S*s for S in SZ for s in sz]
    end
    return sparse(Diagonal(SZ))
end

function generate_string(N::Int64)
    sz = sparse([1 0; 0 -1])
    SZ = 1
    for n=1:N
        SZ = SZ ⊗ sz
    end
    return SZ
end

function generate_fermion_ops(s::Float64, N::Int64)
    spin_matrices = generate_spin_matrices(s)
    d = size(spin_matrices["Id"],1)

    Ô = Dict()

    rn = ["Adag", "A"]
    ρ = real(spin_matrices["S+"]*spin_matrices["S-"])
    N̂ = spzeros(Float64, d^N, d^N)
    sz = 2*spin_matrices["Sz"]
    idL = 1
    for n=1:N
        # @show n
        idR = sparse(I, d^(N-n), d^(N-n))
        for (ids,sm) in enumerate(["S+", "S-"])
            mat = idL ⊗ spin_matrices[sm] ⊗ idR
            Ô[rn[ids], n] = mat
            Ô["Id", n] = sparse(I, d^N, d^N)
        end
        idL = idL ⊗ -sz
        idL2 = sparse(I, d^(n-1), d^(n-1))
        N̂ += idL2 ⊗ ρ ⊗ idR
    end
    return Ô, N̂
end

function are_fermionic_operators(Ô)
    is_fermion = true
    for k1 in keys(Ô), k2 in keys(Ô)
        if k1[1] != "Id" && k2[1] != "Id"
            ac = Ô[k1]*Ô[k2]+Ô[k2]*Ô[k1]
            if length(nonzeros(ac)) != 0
                is_identity = (ac ≈ sparse(I, size(ac)))
                is_fermion = !is_identity ? false : true
            end
        end
    end
    return true
end