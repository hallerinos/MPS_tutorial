using ITensors, BenchmarkTools, LinearAlgebra, KrylovKit
include("bonds.jl")
include("bonds_plotter.jl")
include("lattices.jl")
include("hamiltonians.jl")
include("contract_MPOS.jl")
include("tensor_to_matrix.jl")
include("plot_magnetization.jl")

Nx, Ny = 2, 2
N, graph = square(Nx, Ny)  # see available graphs in lattices.jl
conserve_qns = false
plot_graph2d(graph)

t = 1.0; U=2.0; V=10.0;
sites = siteinds("Electron", N; conserve_qns=conserve_qns)
mpos = Hubbard_2D(graph, sites; t=t, U=U, V=V)

values = nothing
hamiltonian_tensor = nothing
if N < 6  # only compare with ED if less than 8 sites
    @time hamiltonian_tensor = contract_MPOS(mpos)
    @show size(hamiltonian_tensor)
    @time hamiltonian_matrix = tensor_to_matrix(hamiltonian_tensor)
    @time values, vectors = eigsolve(hamiltonian_matrix, 5, :SR, eltype(hamiltonian_matrix), issymmetric=true);
end

# --------------- MPS settings ----------------
n_ex = 1  # how many excitations
M = 64  # set the maximum bond dimension
Ns = 100  # set the maximum number of sweeps
etresh = 1e-12  # naïve stopping criterion
restart = false  # restarting from states MPS
outputlevel = 1  # increase output from 0,1,2
state = [isodd(n) ? "UpDn" : "0" for n in 1:N]  # half filling
psi = MPS(sites, state)  # create random initial state
plot_density_electrons(psi, graph);

sweeps = Sweeps(Ns)  # initialize sweeps object
maxdim!(sweeps, M)  # set bond dimension
obs = DMRGObserver(; energy_tol=etresh)  # observe state after each sweep
psis = Vector{typeof(psi)}(undef, n_ex)  # collect all dmrg states
enes = Vector(undef, n_ex)  # collect all dmrg energies
for n=1:n_ex
    ene, psi = dmrg(mpos, psis[1:n-1], psi, sweeps, weight=1000, observer=obs, outputlevel=outputlevel)
    psis[n] = psi; enes[n] = ene;  # add state and energy to arrays
end
if values != nothing
    @show enes .- values[1:n_ex]  # deviation with respect to ED
end

[plot_density_electrons(psi, graph) for psi in psis];