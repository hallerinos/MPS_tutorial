using BenchmarkTools, KrylovKit
include("../bonds.jl")
include("../bonds_plotter.jl")
include("../lattices.jl")
include("hamiltonians.jl")
include("generate_lattice_ops.jl")

Nx, Ny = 16, 1
N, graph = square(Nx, Ny)  # see available graphs in lattices.jl
# plot_graph2d(graph)

tol = eps()
nev = 5
kdmin = 10  # minimum krylovdim
exact_energy = 1 - 1/sin(π/(4*N+2))  # energy at criticality

@time Ô, N̂ = generate_spin_ops(0.5, N)
@time Ĥ = generate_Ising(Ô, graph)

@time Es, Ψs = eigsolve(Ĥ, nev, :SR, eltype(Ĥ), issymmetric=true, krylovdim=max(nev,kdmin), tol=tol)
@show 2*Es[1] - exact_energy
plot(Es)
;
