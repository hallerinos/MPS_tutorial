include("iDMRG.jl")
include("measure.jl")
using Plots

# loop over different magnetic fields
hs = 0:0.05:2
mz = zeros(length(hs))
en = zeros(length(hs))
for (idh, h) in enumerate(hs)
    ψ, ε = iDMRG(J=1., h=-h, verbose=0, n=20)
    corrs, polarization = measure_obs(ψ)
    mz[idh] = polarization[end]
    en[idh] = ε
end
display(plot(hs, mz, xlabel = "h/J", ylabel="Mz"))