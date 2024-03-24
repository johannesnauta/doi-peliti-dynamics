#= Module for plotting =#
#/ Start module
module Plotter

#/ Packages
using CairoMakie
using FHist
using LaTeXStrings

using Symbolics

#################
### FUNCTIONS ###
function plot_hists(hfull, hred)
    fig = CairoMakie.Figure(; size = (3.25, 3.25/4*3) .* 72)

    edges = collect(hfull.hist.edges[begin])
    limits = (minimum(edges), maximum(edges), 0, 0.02)
    xplot = edges[begin:end-1] .+ diff(edges) / 2
    
    ax = Axis(
        fig[1,1],
        xlabel = L"x_1", ylabel = L"P(x_1)",
        xlabelsize = 15, ylabelsize = 15, xticklabelsize = 10, yticklabelsize = 10,
        xgridvisible=false, ygridvisible=false,        
        xticks = LinRange(limits[1]:50:limits[2]), yticks = LinRange(0, 0.02, 5),
        limits = limits,
    )

    colors = (:black, :rebeccapurple)
    markers = (:circle, :utriangle)
    
    stairs!(
        ax, xplot, hfull.hist.weights, color=:black,
        gap = 0, linewidth = .9,
        label = L"\textrm{full}", label_size=1,
    )
    scatter!(
        ax, xplot, hred.hist.weights,
        color=(:white,0), marker = :diamond, markersize=7,
        strokecolor=:rebeccapurple, strokewidth=1.,
        label = L"\textrm{reduced}"
    )

    #/ Add legend
    axislegend(
        ax, merge=true, position=:lt, patchsize = (5, 5),
        labelsize=12, padding=5.0, framevisible=false, rowgap=0,
        margin = (0.0, 0.0, 0.0, 0.0)
    )
    return fig
end

function get_stationary()
    @variables u
    Du = Differential(u)
    
    f = 1/2 + u/2 -
        ((425 - 3*sqrt(70) + 5*sqrt(7315 - 102*sqrt(70))) * u) / 10800 -
        u^2 / 1000 -
        (u * (500 - u + sqrt(252200 - 1000 * u + u^2))) / 33000 -
        (u * (
            500 - 5/16 * (1000 - u + sqrt(1006400 - 2000 * u + u^2)) +
                sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000 * u + u^2)))^2)
        )) / 20000

    g = 0.5*sqrt(
        1/2 + u/2 +
        ((425 - 3*sqrt(70) + 5*sqrt(7315 - 102*sqrt(70))) * u) / 10800 +
        u^2 / 1000 -
        (u * (u + (1000*u - 2*u^2) / (2*sqrt(252200 - 1000*u + u^2)))) / 16500 +
        (u * (500 - u + sqrt(252200 - 1000*u + u^2))) / 33000 -
        (u * (
            -(5/16) * (u + (2000*u - 2*u^2) / (2*sqrt(1006400 - 2000*u + u^2))) +
                (5 * (u + (2000*u - 2*u^2) / (2*sqrt(1006400 - 2000*u + u^2))) *
                (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))) /
                (16 * sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))^2))
        )) / 10000 +
            (u * (
                500 - 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)) +
                    sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))^2)
            )) / 20000
    )^2

    dg = expand_derivatives(Du(g))
    F_fun = eval(build_function(f, u))
    D_fun = eval(build_function(g, u))
    dD_fun = eval(build_function(dg, u))
    return F_fun, D_fun, dD_fun
end

end # module Plotter
#/ End module
