using ITensors, BenchmarkTools, LinearAlgebra
include("bonds.jl")
include("bonds_plotter.jl")
include("lattices.jl")
include("hamiltonians.jl")
include("contract_MPOS.jl")
include("tensor_to_matrix.jl")
include("plot_local_obs.jl")

Nx, Ny = 3, 2
N, graph = square(Nx, Ny)  # see available graphs in lattices.jl
# plot_graph2d(graph)
conserve_qns = true

t = 1.0; V=1.0;
sites = siteinds("Fermion", N; conserve_qns=conserve_qns)
mpos = Hubbard_spinless(graph, sites; t=t, V=V)

values = nothing
hamiltonian_tensor = nothing
if N < 7  # only compare with ED if less than 8 sites
    @time hamiltonian_tensor = contract_MPOS(mpos)
    @show size(hamiltonian_tensor)
    @time hamiltonian_matrix = tensor_to_matrix(hamiltonian_tensor)
    @time values, vectors = eigen(hamiltonian_matrix)
end

# --------------- MPS settings ----------------
n_ex = 2  # how many excitations
M = 256  # set the maximum bond dimension
Ns = 100  # set the maximum number of sweeps
etresh = 1e-12  # naïve stopping criterion
restart = false  # restarting from states MPS
outputlevel = 1  # increase output from 0,1,2
state = [isodd(n) ? "1" : "0" for n in 1:N]  # half filling
# state[2] = "1"
psi = conserve_qns ? MPS(sites, state) : randomMPS(sites, M)
plot_density_spinless_fermions(psi, graph)

sweeps = Sweeps(Ns)  # initialize sweeps object
maxdim!(sweeps, M)
obs = DMRGObserver(; energy_tol=etresh)
psis = Vector{typeof(psi)}(undef, n_ex)
enes = Vector(undef, n_ex)
for n=1:n_ex
    ene, psi = dmrg(mpos, psis[1:n-1], psi, sweeps, weight=1000, observer=obs, outputlevel=outputlevel)
    psis[n] = psi; enes[n] = ene;  # add state and energy to arrays
end
if values != nothing
    @show enes
    @show enes .- values[1:n_ex]  # deviation with respect to ED
end

[plot_density_spinless_fermions(psi, graph) for psi in psis];