function Ising_MPO(graph::Graph, sites; h::Float64=1.0, J::Float64=1.0)
    # automated MPO generation by a sum of operator expressions
    ampo = OpSum()
    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])

    𝐒 = ["Sx", "Sy", "Sz"]  # vector denoting the spin matrices 0.5σᵢ
    # uniform magnetic field on all sites
    for n in nodes
        ampo += 2*h, 𝐒[3], n[1]
    end
    # loop over all bonds of the graph
    for b in graph
        ampo += 4*J, 𝐒[1], b.s1, 𝐒[1], b.s2
    end
    return MPO(ampo,sites)
end

function Hubbard(graph::Graph, sites; t::Float64=1.0, U::Float64=0.0, V::Float64=0.0)
    # automated MPO generation by a sum of operator expressions
    ampo = OpSum()
    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])

    # on-site Hubbard interaction on all sites
    for n in nodes
        ampo += U, "Nupdn", n[1]
    end
    Cdag = ["Adagup", "Adagdn"]  # vector denoting the spin matrices 0.5σᵢ
    C = ["Aup", "Adn"]  # vector denoting the spin matrices 0.5σᵢ
    # loop over all bonds of the graph
    for b in graph
        for s=1:length(Cdag)
            ampo += -t, Cdag[s], b.s1, C[s], b.s2
            ampo += -t, Cdag[s], b.s2, C[s], b.s1
        end
        ampo += V, "Ntot", b.s1, "Ntot", b.s2
    end
    return MPO(ampo,sites)
end

function Hubbard_spinless(graph::Graph, sites; t::Float64=1.0, V::Float64=0.0)
    # automated MPO generation by a sum of operator expressions
    ampo = OpSum()
    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])

    # loop over all bonds of the graph
    for b in graph
        ampo += -t, "Cdag", b.s1, "C", b.s2
        ampo += -t, "Cdag", b.s2, "C", b.s1
        ampo += V, "N", b.s1, "N", b.s2
    end
    return MPO(ampo,sites)
end