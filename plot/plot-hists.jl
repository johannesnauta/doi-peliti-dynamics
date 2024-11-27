#= Module for plotting histograms and related plots =#
#/ Start module
module HistPlotter

#/ Packages
using CairoMakie
using FHist
using JLD2
using LaTeXStrings

#################
### FUNCTIONS ###
"""
 Plot histogram of reduced (stochastic) dynamics
 If desired, overlays histogram of full (stochastic) dynamics (if it exists)
"""
function plot_rhists(
    ;
    hfname = "../data/pdfs/pdf_rSIRS.jld2",
    width = .7*246,
    overlay = true,
    savefig = false,
    figname = nothing 
    )
    #/ Create figure
    fig = Figure(
        ; size=(width,width/1.5), figure_padding=(2,2,2,4), backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1],
        xlabel=L"n_S", ylabel=L"p(n_S)", xlabelsize=12, ylabelsize=12,
        xticklabelsize=8, yticklabelsize=8,
        xlabelpadding=0.0,
        yticks=[0,2,4,6],
        xminorticksvisible=true, xminorticks=IntervalsBetween(5),
        yminorticksvisible=true, yminorticks=IntervalsBetween(2),
        xgridvisible=false, ygridvisible=false,
        limits=(3.0, 4.5, 0.0, 7.0)
    )

    #/ Load data
    rdb = JLD2.load(hfname)
    rfh = rdb["histogram"]
    #/ Plot
    edges = rfh.binedges[begin].nonuniform_edges
    barpositions = (edges[1:end-1] .+ edges[2:end]) ./ 2
    br = barplot!(
        barpositions ./ 1e3, rfh.bincounts * 1e3, color=:rebeccapurple, gap=0,
        label=L"\textrm{reduced}"
    )

    #/ Overlay histogram of full dynamics as stair plot
    if overlay
        #~ Extract filename from given filename of reduced histogram
        fname = replace(hfname, "_r" => "_")
        db = JLD2.load(fname)
        fh = db["histogram"]
        _edges = rfh.binedges[begin].nonuniform_edges
        stairpositions = (edges[1:end-1] .+ edges[2:end]) ./ 2
        sr = stairs!(
            stairpositions ./ 1e3, fh.bincounts * 1e3, color=:black, step=:center,
            linewidth=0.8, label=L"\textrm{full}"
        )
        axislegend(
            ax, position=:rt, framevisible=false, labelsize=8.5, patchsize=(4,4),
            padding=0
        )
    end

    #/ Some labels
    Label(fig[1,1,Top()], halign=:left, L"\times 10^{-3}", fontsize=8)
    Label(fig[1,1,Right()], valign=:bottom, halign=:right, L"\times 10^{3}", fontsize=8)

    return fig
end

end # module HistPlotter
#/ End module
