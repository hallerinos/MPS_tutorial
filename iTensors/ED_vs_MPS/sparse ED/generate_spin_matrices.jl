using SparseArrays
function generate_spin_matrices(physicalSpin::Float64)

    # initialize spin matrices
    d = Int(2 * physicalSpin + 1);
    Sx = zeros(ComplexF64,d,d);
    Sy = zeros(ComplexF64,d,d);
    Sz = zeros(ComplexF64,d,d);

    # construct Sx, Sy and Sz
    for idx = 1 : d

        for idy = 1 : d

            entryXY = 0.5 * sqrt((physicalSpin + 1) * (idx + idy - 1) - idx * idy);

            if (idx + 1) == idy
                Sx[idx,idy] += entryXY;
                Sy[idx,idy] -= 1im * entryXY;
            end

            if idx == (idy + 1)
                Sx[idx,idy] += entryXY;
                Sy[idx,idy] += 1im * entryXY;
            end

            if idx == idy
                Sz[idx,idy] += physicalSpin + 1 - idx;
            end

        end

    end

    # compute Id, Sp and Sm
    Id = one(Sz);
    Sp = Sx + 1im * Sy;
    Sm = Sx - 1im * Sy;

    spin_matrices = Dict("Sx" => sparse(Sx), "Sy" => sparse(Sy), "Sz" => sparse(Sz), "S+" => sparse(Sp), "S-" => sparse(Sm), "Id" => sparse(Id))

    # return spin matrices
    return spin_matrices

end