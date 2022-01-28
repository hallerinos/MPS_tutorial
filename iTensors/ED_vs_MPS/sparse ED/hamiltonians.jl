function generate_Ising(Ô::Dict, graph::Graph; h::Float64=1.0, J::Float64=2.0)
    Ĥ = spzeros(ComplexF64, size(Ô["Id",1])...)

    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])

    # uniform magnetic field
    for n in nodes
        Ĥ .+= h*Ô["Sz", n[1]]
    end

    # nearest neighbor interactions
    for b in graph
        Ĥ .+= J*Ô["Sx", b.s1]*Ô["Sx", b.s2]
    end
    return Ĥ
end

function generate_Hubbard(Ô::Dict, graph::Graph; t::Float64=1.0, V::Float64=0.0, μ::Float64=0.0)
    Ĥ = spzeros(Float64, size(Ô["Id",1])...)

    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    for n in nodes
        # @show n
        Ĥ .+= -μ.*(Ô["Adag", n[1]]*Ô["A", n[1]])
    end
    # nearest neighbor interactions
    for b in graph
        # @show b
        Ĥ .+= -t.*(Ô["Adag", b.s1]*Ô["A", b.s2])
        Ĥ .+= -t.*(Ô["Adag", b.s2]*Ô["A", b.s1])

        Ĥ .+= V.*((Ô["Adag", b.s1]*Ô["A", b.s1])*(Ô["Adag", b.s2]*Ô["A", b.s2]))
    end
    return Ĥ
end