#= Module for plotting trajectories =#
#/ Start module
module TPlotter

#/ Packages
using CairoMakie
using JLD2
using LaTeXStrings

#################
### FUNCTIONS ###
function plot_trajectories_mSIRS(;
    ddir = "../data/trajectories/",
    fname = ddir*"mSIRStrajectories.jld2",
    width = .7*246,
    norm = 1e4,
    savefig = false,
    figname = nothing
)
    #/ Create figure
    fig = Figure(
        ; size=(width, width/1.77), figure_padding=(2,2,2,4), backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1], xlabel=L"t", ylabel=L"n_S", xlabelsize=12, ylabelsize=12,
        xticklabelsize=8, yticklabelsize=8,
        xlabelpadding=0.0, yticks=[0.0, 1.0],
        xminorticksvisible=true, xminorticks=IntervalsBetween(2),
        yminorticksvisible=true, yminorticks=IntervalsBetween(5),
        xgridvisible=false, ygridvisible=false,
        limits=(0.0, 5.0, 0.0, 1.0)
    )

    #/ Load data
    db = JLD2.load(fname)
    xsde = db["xsde"]
    tsde = db["tsde"]
    xode = db["xode"]
    tode = db["tode"]
    xr   = db["xr"]
    tr   = db["tr"]

    #/ Plot
    xv = [xsde, xr, xode]
    tv = [tsde, tr, tode]
    linewidths = [0.65, 0.65, 0.8]
    colors = [:black, :rebeccapurple, :firebrick2]
    linestyles = [:solid, :solid, (:dash,:dense)]
    labels = [L"\textrm{full}", L"\textrm{reduced}", L"\textrm{mean‐field}"]
    for i in eachindex(xv)
        lines!(
            ax, tv[i] ./ 1e1, xv[i] ./ norm, color=colors[i], label=labels[i],
            linestyle=linestyles[i], linewidth=linewidths[i]
        )
    end
    axislegend(
        ax, position=:rt, labelsize=8.5, patchsize=(6,1), patchlabelgap=2,
        rowgap=0, colgap=0, nbanks=1, 
        padding=0, margin=(0,2,0,2), framevisible=false
        
    )

    #/ Some labels
    Label(fig[1,1,Top()], halign=:left, L"\times 10^{4}", fontsize=8)
    Label(fig[1,1,Right()], valign=:bottom, halign=:right, L"\times 10^{1}", fontsize=8)

    (savefig && !isnothing(figname)) && (save(figname, fig, pt_per_unit=1))
    return fig
end

function plot_trajectories_gLV(;
    ddir = "../data/trajectories/",
    fname = ddir*"gLVtrajectories.jld2",
    width = .7*246,
    norm = 1e3,
    savefig = false,
    figname = nothing
)
    #/ Create figure
    fig = Figure(
        ; size=(width, width/1.77), figure_padding=(2,2,2,4), backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1], xlabel=L"t", ylabel=L"n_1", xlabelsize=12, ylabelsize=12,
        xticklabelsize=8, yticklabelsize=8,
        xlabelpadding=0.0, yticks=[0.0, 1.0],
        xminorticksvisible=true, xminorticks=IntervalsBetween(2),
        yminorticksvisible=true, yminorticks=IntervalsBetween(5),
        xgridvisible=false, ygridvisible=false,
        limits=(0.0, 5.0, 0.0, 1.05)
    )

    #/ Load data
    db = JLD2.load(fname)
    xsde = db["xsde"]
    tsde = db["tsde"]
    xode = db["xode"]
    tode = db["tode"]
    xr   = db["xr"]
    tr   = db["tr"]

    #/ Plot
    xv = [xsde, xr, xode]
    tv = [tsde, tr, tode]
    linewidths = [0.65, 0.65, 0.8]
    colors = [:black, :rebeccapurple, :firebrick2]
    linestyles = [:solid, :solid, (:dash,:dense)]
    labels = [L"\textrm{full}", L"\textrm{reduced}", L"\textrm{mean‐field}"]
    for i in eachindex(xv)
        lines!(
            ax, tv[i] ./ 1e1, xv[i] ./ norm, color=colors[i], label=labels[i],
            linestyle=linestyles[i], linewidth=linewidths[i]
        )
    end
    axislegend(
        ax, position=:rb, labelsize=8.5, patchsize=(6,1), patchlabelgap=2,
        rowgap=0, colgap=0, nbanks=1, 
        padding=0, margin=(0,2,2,0), framevisible=false
        
    )

    #/ Some labels
    Label(fig[1,1,Top()], halign=:left, L"\times 10^{3}", fontsize=8)
    Label(fig[1,1,Right()], valign=:bottom, halign=:right, L"\times 10^{1}", fontsize=8)

    (savefig && !isnothing(figname)) && (save(figname, fig, pt_per_unit=1))
    return fig
end


function plot_trajectories_sLV(;
    Dv=[0.1, 0.5, 0.9],
    ddir="../data/trajectories/",    
    savefig = false,
    figname = nothing
)
    #/ Define figure and axes
    width = .8*246
    fig = Figure(; size=(width, width), figure_padding=(2,8,2,4), backgroundcolor=:transparent)
    axes = [Axis(
        fig[i,1], xlabel=L"t", ylabel=L"\langle n_i \rangle",
        xlabelsize=12, ylabelsize=12,
        xlabelvisible=(i==3), ylabelvisible=true,
        xgridvisible=false, ygridvisible=false,
        xticklabelsize=8, yticklabelsize=8,
        yticks = [0.0, 0.5],
        xticklabelsvisible=(i==3),
        xticksize=3, xminorticksize=1.5, xminorticksvisible=true,
        yticksize=3, yminorticksize=1.5, yminorticksvisible=true,
        limits = (0,50,0.0,0.6)
    ) for i in eachindex(Dv)]

    colors = [:black, :rebeccapurple]
    labels = [L"\langle n_A \rangle", L"\langle n_B \rangle"]
    dtext = [L"D<D_1", L"D_1 < D < D_2", L"D > D_2"]

    for i in eachindex(Dv)
        #/ Load data
        fname = ddir * "sLVtrajectory_D$(Dv[i]).jld2"
        db = JLD2.load(fname)
        x = db["x"]
        t = db["t"]
        #/ Plot
        for j in 1:2
            l = lines!(
                axes[i], t, x[j,:], color=colors[j], linewidth=.8, label=labels[j]
            )
        end
        txt = text!(
            axes[i], 0.95, 0.95, text=dtext[i], space=:relative, align=(:right,:top),
            fontsize=10
        )
    end

    #/ Legend
    axislegend(
        axes[end], framevisible=false, labelsize=10, position=:lt, padding=0,
        patchsize=(5,1), rowgap=0, margin=(4,0,0,1), patchlabelgap=1
    )

    #/ Adjust
    rowgap!(fig.layout, 8)
    resize_to_layout!(fig)

    (savefig && !isnothing(figname)) && (save(figname, fig, pt_per_unit=1))
    return fig
end


end # module Plotter
#/ End module
