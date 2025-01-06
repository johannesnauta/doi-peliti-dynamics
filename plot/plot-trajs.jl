#= Module for plotting trajectories =#
#/ Start module
module TPlotter

#/ Packages
using CairoMakie
using JLD2
using LaTeXStrings

#################
### FUNCTIONS ###
function plot_trajectories_sLV(;
    Dv=[0.1, 0.5, 0.9],
    ddir="../data/trajectories/"
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
    
    return fig
end


end # module Plotter
#/ End module
