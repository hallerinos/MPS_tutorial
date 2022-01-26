# remove matplotlib dependency here...
using PyPlot, PyCall
pygui(true)

function plot_graph(graph::Graph)
    len = unique([length(b.r1) for b in graph])
    check = any([length(b.r2) == length(b.r1) for b in graph])
    @show length(len), check
    if !check || length(len)>1
        throw(ArgumentError("r1 and r2 have incompatible shapes"))
    end
    if len[1]<=2
        return plot_graph2d(graph)
    elseif len[1]==3
        return plot_graph3d(graph)
    else
        throw(ArgumentError("no plotting function availavble for lattice dimensions > 3"))
    end
end

"""
Plot graph in 3d
"""
function plot_graph3d(graph::Graph)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    x = [b.r1[1] for b in graph]
    y = [b.r1[2] for b in graph]
    z = [b.r1[3] for b in graph]
    u = [b.r2[1] - b.r1[1] for b in graph]
    v = [b.r2[2] - b.r1[2] for b in graph]
    w = [b.r2[3] - b.r1[3] for b in graph]
    ax[:quiver](x,y,z, u,v,w)
    xlabel(L"x")
    ylabel(L"y")
    zlabel(L"z")
    return fig
end

"""
Plot x,y plane structure
"""
function plot_graph2d(graph::Graph)
    # plt.rc("text", usetex=true)
    x = [b.r1[1] for b in graph]
    y = [b.r1[2] for b in graph]
    u = [b.r2[1] - b.r1[1] for b in graph]
    v = [b.r2[2] - b.r1[2] for b in graph]
    fig = plt.quiver(x,y,u,v, angles="xy", scale_units="xy", scale=1)
    ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    plt.scatter([n[2][1] for n in ns],[n[2][2] for n in ns])
    xlabel(L"x")
    ylabel(L"y")
    return fig
end