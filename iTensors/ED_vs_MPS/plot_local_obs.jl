function plot_magnetization(psi, graph)
    ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    xs = [n[2][1] for n in ns]
    ys = [n[2][2] for n in ns]
    lobs = [expect(psi, lob) for lob in ["Sx", "Sy", "Sz"]]
    fig = figure(figsize=2.0.*((3+3/8), (3+3/8)/1.25))
    scatter(xs, ys, cmap="RdBu_r", c=lobs[3], marker="h", s=1600, vmin=-0.5, vmax=0.5)
    quiver(xs, ys, lobs[1],lobs[2], scale=.5, units="xy", pivot="middle", color="white")
    maxx = maximum(xs)
    maxy = maximum(ys)
    minx = minimum(xs)
    miny = minimum(ys)
    dy = (maxy-miny)/Ny
    dx = (maxx-minx)/Nx
    ylim(miny-dy,maxy+dy)
    xlim(minx-dx,maxx+dx)
    axis("off")
    return fig
end

function plot_density_electrons(psi, graph)
    ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    xs = [n[2][1] for n in ns]
    ys = [n[2][2] for n in ns]
    lobs = [expect(psi, lob) for lob in ["Ntot", "Nup", "Ndn"]]
    fig = figure(figsize=2.0.*((3+3/8), (3+3/8)/1.25))
    scatter(xs, ys, c=lobs[1], marker="h", s=1600, vmin=0, vmax=2)
    quiver(xs, ys, 0.0*lobs[2], lobs[2]-lobs[3], units="xy", pivot="middle", color="white")
    maxx = maximum(xs)
    maxy = maximum(ys)
    minx = minimum(xs)
    miny = minimum(ys)
    dy = (maxy-miny)/Ny
    dx = (maxx-minx)/Nx
    ylim(miny-dy,maxy+dy)
    xlim(minx-dx,maxx+dx)
    axis("off")
    return fig
end