function contract_MPOS(mpos::MPO)::ITensor
    hamiltonian_tensor = mpos[1]
    for mpo in mpos[2:end]
        @disable_warn_order hamiltonian_tensor *= mpo
    end
    return hamiltonian_tensor
end