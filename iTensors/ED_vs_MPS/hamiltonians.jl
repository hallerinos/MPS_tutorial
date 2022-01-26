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
        dir = b.r2 .- b.r1  # the direction vector between lattice nodes
        dist = norm(dir)
        if dist ≈ 1  # nearest neighbor interaction
            ampo .+= 4*J, 𝐒[1], b.s1, 𝐒[1], b.s2
        end
    end
    return MPO(ampo,sites)
end