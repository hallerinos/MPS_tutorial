using BenchmarkTools, KrylovKit
include("../bonds.jl")
include("../bonds_plotter.jl")
include("../lattices.jl")
include("hamiltonians.jl")
include("generate_lattice_ops.jl")

Nx, Ny = 3, 2
N, graph = square(Nx, Ny; xperiodic=false, yperiodic=false)  # see available graphs in lattices.jl
# plot_graph2d(graph)

tol = eps()
nev = 7
kdmin = 10  # minimum krylovdim
n = N/2-1

@time Ô, N̂ = generate_fermion_ops(0.5, N)
pind = getindex.(findall(x -> x == n, N̂), 1)

# @time @show are_fermionic_operators(Ô)
@time Ĥ = generate_Hubbard(Ô, graph; t=-1.0, μ=0.0, V=1.0)
@time Ĥₙ = Ĥ[pind,pind]

@time Es, Ψs = eigsolve(Ĥₙ, nev, :SR, eltype(Ĥ), issymmetric=true, krylovdim=max(nev,kdmin), tol=tol);
@show Es