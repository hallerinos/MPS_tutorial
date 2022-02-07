# remove matplotlib dependency here...
using PyCall
@pyimport matplotlib.pyplot as plt
@pyimport numpy as np
plt.rc("text", usetex=true)
plt.rc("font", family="serif")

# just an example
function testplot()
    plt.close()
    fig = plt.figure()
    x = np.linspace(0, 2*np.pi, 1000)
    y = np.sin(3*x + 4 * np.cos(2*x))
    plt.plot(x, y, color="red", linewidth=2.0, linestyle="--")
    plt.title("A sinusoidally modulated sinusoid")
    # plt.savefig("test.pdf")
    return fig
end

function plot_graph(graph::Graph)
    plt.close()
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
    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    x = [b.r1[1] for b in graph]
    y = [b.r1[2] for b in graph]
    z = [b.r1[3] for b in graph]
    u = [b.r2[1] - b.r1[1] for b in graph]
    v = [b.r2[2] - b.r1[2] for b in graph]
    w = [b.r2[3] - b.r1[3] for b in graph]
    ax[:quiver](x,y,z, u,v,w)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.zlabel("z")
    return fig
end

"""
Plot x,y plane structure
"""
function plot_graph2d(graph::Graph)
    plt.close()
    x = [b.r1[1] for b in graph]
    y = [b.r1[2] for b in graph]
    u = [b.r2[1] - b.r1[1] for b in graph]
    v = [b.r2[2] - b.r1[2] for b in graph]
    fig = plt.quiver(x,y,u,v, angles="xy", scale_units="xy", scale=1)
    ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    plt.scatter([n[2][1] for n in ns],[n[2][2] for n in ns])
    plt.xlabel("x")
    plt.ylabel("y")
    return fig
end