#/ Start module
module PhasePlotter

#/ Packages
using CairoMakie
using LaTeXStrings

#################
### FUNCTIONS ###
"Plot η-value (occupation no.) of global fixed point F vs γ, for both α"
function plot_occupation_mSIS(; α=0.5, ρ=1.0, N=1.0)
    #/ Create figures
    width = .8*246
    fig = Figure(
        ; size=(width,.6*width), figure_padding=(1,6,1,3),
        backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1],
        xlabel=L"\textrm{control\;parameter}\;\gamma",
        ylabel=L"\textrm{order\;parameter}\;\varphi",
        xlabelsize=11, ylabelsize=11,
        xlabelpadding=1.0, ylabelpadding=3.0,
        yticks=[0.0, 1.0], xticks=[0.0, 1.0, 2.0],
        xticklabelsize=8, yticklabelsize=8,
        xminorticksvisible=true, xminorticks=IntervalsBetween(2),
        yminorticksvisible=true, yminorticks=IntervalsBetween(2),
        limits=(0.0, 2ρ, 0., 1.01),
        xgridvisible=false, ygridvisible=false
    )
    #/ Define function for the order parameter (from theory, see below)
    ϕ(γ,α) = 1 - ((α+γ)*N - sqrt((α+γ)^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)
    #~ specify γ values
    γv = ρ:0.01:2ρ
    αlabel = α < ρ ? L"\alpha < \rho" : L"\alpha > \rho"
    color = α < ρ ? :firebrick2 : :firebrick2

    #/ Plot lines
    lines!(ax, vcat(0.0, ρ), vcat(0.0, 0.0), color=color, linewidth=1.2, label=αlabel)
    lines!(ax, γv, ϕ.(γv,α), color=color, linewidth=1.2)
    #~ if 1st order transition, identify it more clearly    
    if α > ρ
        vlines!(
            ax, [1.0], ymin=0.0, ymax=ϕ(γv[begin], α), color=:gray,
            linewidth=1., linestyle=:dot
        )
    end

    axislegend(
        ax, position=:rb, labelsize=11, framevisible=false, rowgap=0,
        patchsize=(0,0), padding=0
    )
    return fig 
end

"Plot some arrowheads at all (x,y) with the direction given by the field"
function plot_arrowfield(ax, field, x, y, color)
    derivatives = field.(x, y)
    dx, dy = first.(derivatives), last.(derivatives)
    # arrowpoly = Makie.Polygon(Point2f[(0, 1), (0.5, -0.3), (0, 0), (-0.5, -0.3)])
    arrows!(
        ax, x, y, dx, dy, color=color,
        linewidth=0., align=:center, arrowsize=7
    )
    if color == :rebeccapurple
        arrows!(
            ax, x, y, dx, dy, color=:white,
            linewidth=0., align=:center, arrowsize=3.5
        )
    end
end

function plot_phase_SIR(;
    N::Float64=1.0, ρ::Float64=1.0, α::Float64=2.0, γ=[2.0/3.0,1.0,2.0],
    multiple_figures=false
    )
    #~ field copied from Symbolics output
    field(η,θ,γ) = (
        ((N^3)*exp(θ)*ρ - (N^2)*exp(θ)*η*ρ - (N^2)*exp(-θ)*α*η - (N^2)*exp(-θ)*γ*η
         + 2.0N*exp(-θ)*α*(η^2) + N*exp(-θ)*γ*(η^2) - exp(-θ)*α*(η^3)) / (N^2),
        (-(N^2)*ρ + (N^2)*exp(θ)*ρ - (N^2)*exp(-θ)*α - (N^2)*exp(-θ)*γ +
         4.0N*exp(-θ)*α*η + 2.0N*exp(-θ)*γ*η - 3.0exp(-θ)*α*(η^2) +
         (N^2)*exp(θ)*exp(-θ)*α + (N^2)*exp(θ)*exp(-θ)*γ - 4.0N*exp(θ)*exp(-θ)*α*η -
         2.0N*exp(θ)*exp(-θ)*γ*η + 3.0exp(θ)*exp(-θ)*α*(η^2)) / (N^2)
    )
    #~ field copied from doing the derivatives by hand (and/or using Mathematica)
    # field(η,θ,γ) = (
    #     -exp(-θ)*(N-η) * (-N*(α+γ)*η + α*η^2 + exp(2*θ)*N^2*ρ) / N^2,
    #     -exp(-θ)*(-1+exp(θ)) * (-2*N*(2*α + γ)*η + 3*α*η^2 + N^2 * (α+γ + exp(θ)*ρ)) / N^2
    # )
    f(η,θ,γ) = Point2f(field(η,θ,γ))
    
    width = .8*246
    fig = !(multiple_figures) ? Figure(;
        size=(4*width/3,width/1.8), figure_padding=(0,4,2,4),
        backgroundcolor=:transparent
    ) : [Figure(;
        size=(0.33*4*width/3,width/2), figure_padding=(0,0,0,0),
        backgroundcolor=:transparent
    ) for _ in 1:3]
    axtitles = [L"\gamma < \rho", L"\gamma = \rho", L"\gamma > \rho"]
    
    ax = [Axis(
        !(multiple_figures) ? fig[1,i] : fig[i][1,1],
        title=axtitles[i], titlesize=11,
        xlabel=L"\eta", xlabelpadding=0.0,
        ylabel=L"\theta", ylabelpadding=0.0,
        aspect=1,
        xlabelsize=12, ylabelsize=12, ylabelvisible=(i==1),
        xticks = [0.0, 1.0, 2.0], xticksize=2,
        xminorticksvisible=true, xminorticks=IntervalsBetween(2), xminorticksize=1,
        yticks=[-0.5, 0.0, 0.5, 1.0], yticklabelsvisible=(i==1), yticksize=2,
        xticklabelsize=7, yticklabelsize=7,
        xgridvisible=false, ygridvisible=false,
        limits=(0,2,-.5,1),
    ) for i in 1:length(γ)]
    
    # αlabel = (α > ρ) ? L"\alpha > \rho" : L"\alpha < \rho"
    # γlabel = [L"\gamma<\rho", L"\gamma=\rho", L"\gamma>\rho"]
    
    for i in eachindex(γ)    
        band!(ax[i], [1,2], [-0.5,-0.5], [1.0,1.0], color=(:red, 0.2))
        #/ Add streamplot
        _f(η,θ) = f(η,θ,γ[i])
        sp = streamplot!(
            ax[i], _f,
            0.0..2.0, -0.5..1.0,
            color= x -> :gray, alpha=0.8,
            arrow_size=3.5, linewidth=.3,
            density=.33, maxsteps=500, stepsize=0.005,
            gridsize=(32,32)
        )
        # / Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:black, linewidth=1.2)
        hlines!(ax[i], [0.0], color=:black, linewidth=1.2)
        
        #/ Add non-trivial zero-energy line
        ηmax = min(2.0, N*(α+γ[i])/α - 1e-7)
        ηplot = 0.0:0.01:ηmax
        θf(η) = @. log(η * (α*(N-η) + γ[i]*N)/(ρ*N^2))
        lp = lines!(
            ax[i], ηplot, θf(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )

        #/ Add arrows
        #~ field on trivial zero-energy lines with θ=0
        ηt = γ[i] < ρ ? [0.5, 1.5] : ((γ[i] == ρ) ? [0.75, 1.5] : [0.7, 1.35])
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ηt, zeros(length(ηt)), :black)
        θt = γ[i] < ρ ? [-0.2, 0.5] : ((γ[i] == ρ) ? [-0.3, 0.5] : [-0.3, 0.35])
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ones(length(θt)), θt, :black)
        #~ field on non-trivial zero-energy line(s)
        if γ[i] < ρ
            ηnt = [0.65]
        elseif γ[i] == ρ
            ηnt = [0.3, 0.75, (1.0+ηmax)/2.1]
        else
            ηnt = [0.3, 0.65, (1.0+ηmax)/2.2]
        end
        # ηnt = γ[i] ≤ ρ ? [0.4, 0.75, (1.0+ηmax)/2.1] : [0.2, 0.6, (1.0+ηmax)/2]
        θnt = θf.(ηnt)
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ηnt, θnt, :rebeccapurple)
        
        #/ Scatter fixed points        
        if i == 1
            #~ fixed point F below critical transition
            F = [[1.0, 1.0], [0.0, θf(1.0)]]
            markers = [:star5, :circle]
            markersize = [9, 6]
            colors = [:white, :mediumpurple1]
            strokecolors = [:firebrick2, :rebeccapurple]
        elseif i == 2
            # #~ fixed points at critical transition
            ηmin = ((α+γ[i])*N - sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)
            F = [[1.0, ηmin], [0.0, 0.0]]
            markers = [:circle, :star5]
            markersize = [6,9]
            colors = [:mediumpurple1, :white]
            strokecolors = [:rebeccapurple, :firebrick2]
        else        
            # #~ fixed points above critical transition
            ηmin = ((α+γ[i])*N - sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)
            ηplus = ((α+γ[i])*N + sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)            
            F = [[1.0, ηmin, ηplus, 1.0], [0.0, 0.0, 0.0, θf(1.0)]]
            markers = [:circle, :star5, :circle, :circle]
            markersize = [6,9,6,6]
            colors = [:mediumpurple1, :white, :mediumpurple1, :mediumpurple1]
            strokecolors = [:rebeccapurple, :firebrick2, :rebeccapurple, :rebeccapurple]
        end
        s = scatter!(
            ax[i], F[begin], F[end], marker=markers, markersize=markersize,
            color=colors, strokewidth=.9, strokecolor=strokecolors
        )

        # tα = text!(
        #     ax[i], 0.03, 0.99, text=αlabel, fontsize=12,
        #     space=:relative, align=(:left,:top)
        # )
        # tγ = text!(
        #     ax[i], 0.97, 0.99, text=γlabel[i], fontsize=12,
        #     space=:relative, align=(:right,:top)
        # )
    end

    # gl = GridLayout(fig[1,4], tellwidth=false, tellheight=false, halign=:left)
    # b = Box(gl[1,1], strokecolor=:black, strokewidth=.8, color=:white, cornerradius=5)
    # label = (α > ρ) ? L"\alpha > \rho" : L"\alpha < \rho"
    # lbl = Label(gl[1,1], label, rotation=0, fontsize=12, padding=(4,4,4,4))
    
    # colsize!(fig.layout, 4, Fixed(24))
    colgap!(fig.layout, 0)
    resize_to_layout!(fig)
    return fig
end

function plot_phase_SIS(; N::Float64=1.0, ρ::Float64=1.0, γ = [2.0/3.0,1.0,2.0])
    f(η,θ,γ) = Point2f(
        (N-η) * (γ*η*exp(-θ)/N - ρ*exp(θ)),
        (1 - exp(θ)) * (γ*exp(-θ) - 2*γ*η/N*exp(-θ) + ρ)
    )

    #/ Initialize figure
    width = .8*246
    fig = Figure(;
        size=(4*width/3,width/1.75), figure_padding=(0,0,1,1),
        backgroundcolor=:transparent
    )
    axtitles = [L"\gamma < \rho", L"\gamma = \rho", L"\gamma > \rho"]

    ax = [Axis(
        fig[1,i],
        title=axtitles[i], titlesize=10,
        xlabel=L"\eta", xlabelpadding=0.0,
        ylabel=L"\theta", ylabelpadding=0.0,
        aspect=1,
        xlabelsize=12, ylabelsize=12, ylabelvisible=(i==1),
        xticks = [0.0, 1.0, 2.0], xticksize=2,
        xminorticksvisible=true, xminorticks=IntervalsBetween(2), xminorticksize=1,
        yticks=[-0.5, 0.0, 0.5, 1.0], yticklabelsvisible=(i==1), yticksize=2,
        xticklabelsize=7, yticklabelsize=7,
        xgridvisible=false, ygridvisible=false,
        limits=(0,2,-.5,1),
    ) for i in 1:length(γ)]
    # ax = [Axis(
    #     fig[1,i], title=axtitles[i], titlesize=11,
    #     xlabel=L"\eta",
    #     ylabel=L"\theta", ylabelvisible=(i==1),
    #     xlabelsize=12, ylabelsize=12,
    #     xticklabelsize=7, yticklabelsize=7, yticklabelsvisible=(i==1),
    #     xgridvisible=false, ygridvisible=false,
    #     limits=(0,2,-0.5,1)
    # ) for i in 1:length(γ)]

    strokecolors = [
        [:firebrick2,:black], [:firebrick2,:firebrick2], [:black,:firebrick2]
    ]
    
    for i in eachindex(γ)
        band!(ax[i], [1,2], [-0.5,-0.5], [1.0,1.0], color=(:red, 0.2))
        #/ Add streamplot
        _f(η,θ) = f(η,θ,γ[i])
        sp = streamplot!(
            ax[i], _f,
            0.0..2.0, -0.5..1.0, color= x -> :gray, alpha=.8,
            arrow_size=4., linewidth=.25,
            density=.5, maxsteps=500, stepsize=0.001,
            gridsize=(32,32)
        )
        # / Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:black, linewidth=1.2) 
        hlines!(ax[i], [0.0], color=:black, linewidth=1.2)
        
        #/ Add non-trivial zero-energy line        
        ηplot = 0.01:0.01:2.0
        θf(η) = @. log(η * γ[i] / (ρ * N))
        lp = lines!(
            ax[i], ηplot, θf(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )

        #/ Add arrows
        #~ field on trivial zero-energy lines with θ=-
        ηt = γ[i] == ρ ? [0.5, 1.5] : ((γ[i] < ρ) ? [0.5, 1.25] : [0.75, 1.5])
        plot_arrowfield(ax[i], _f, ηt, zeros(length(ηt)), :black)
        θt = γ[i] ≤ ρ ? [-0.25, 0.5] : [-0.25, 0.25]
        plot_arrowfield(ax[i], _f, ones(length(θt)), θt, :black)
        ηnt = γ[i] < ρ ? [1.25, 1.75] : [0.7, 1.3]
        plot_arrowfield(ax[i], _f, ηnt, θf.(ηnt), :rebeccapurple)
        
        #/ Scatter fixed points        
        if i == 1
            #~ fixed point F below critical transition
            F = [[1.0, 1.0, ρ*N/γ[i]], [0.0, θf(1.0), 0.0]]
            markers = [:star5, :circle, :circle]
            markersize = [9, 6, 6]
            colors = [:white, :mediumpurple1, :mediumpurple]
            strokecolors = [:firebrick2, :rebeccapurple, :rebeccapurple]
        elseif i == 2
            F = [[1.0], [0.0]]
            markers = [:star5]
            markersize = [9]
            colors = [:white]
            strokecolors = [:firebrick2]
        else
            F = [[ρ*N/γ[i], 1.0, 1.0], [0.0, 0.0, θf(1.0)]]
            markers = [:star5, :circle, :circle]
            markersize = [9, 6, 6]
            colors = [:white, :mediumpurple1, :mediumpurple]
            strokecolors = [:firebrick2, :rebeccapurple, :rebeccapurple]
        end
            
            
        s = scatter!(
            ax[i], F[begin], F[end], marker=markers, markersize=markersize,
            color=colors, strokewidth=.9, strokecolor=strokecolors
        )
        # elseif i == 2
        #     # #~ fixed points at critical transition
        #     ηmin = ((γ[i]*N) / (2*α*N)
        #     F = [[1.0, ηmin], [0.0, 0.0]]
        #     markers = [:circle, :star5]
        #     markersize = [6,9]
        #     colors = [:mediumpurple1, :white]
        #     strokecolors = [:rebeccapurple, :firebrick2]
        # else        
        #     # #~ fixed points above critical transition
        #     ηmin = ((α+γ[i])*N - sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)
        #     ηplus = ((α+γ[i])*N + sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)            
        #     F = [[1.0, ηmin, ηplus, 1.0], [0.0, 0.0, 0.0, θf(1.0)]]
        #     markers = [:circle, :star5, :circle, :circle]
        #     markersize = [6,9,6,6]
        #     colors = [:mediumpurple1, :white, :mediumpurple1, :mediumpurple1]
        #     strokecolors = [:rebeccapurple, :firebrick2, :rebeccapurple, :rebeccapurple]
        # η_nontrivial = [1.25, nothing, 0.75]
        # if i!=2
        #     dη, dθ = _f(η_nontrivial[i], θf(η_nontrivial[i]))
        #     arrθ = arrows!(
        #         ax[i], [η_nontrivial[i]], [θf(η_nontrivial[i])], [dη], [dθ],
        #         color=:rebeccapurple, linewidth=0., align=:center, arrowsize=7
        #     )
        #     snontrivial = scatter!(
        #         ax[i], [1.0], [log(γ[i]/ρ)], markersize=6,
        #         color=:white, strokewidth=.9, strokecolor=:rebeccapurple
        #     )
        #     s = scatter!(
        #         ax[i], [1,ρ*N/γ[i]], [0,0], markersize=6,
        #         color=:white, strokewidth=.9, strokecolor=[:black,:firebrick2]
        #     )
        # else
        #     s = scatter!(
        #         ax[i], [1], [0], markersize=6,
        #         color=:white, strokewidth=.9, strokecolor=:firebrick2
        #     )
        # end
        
    end
    colgap!(fig.layout, 0)
    resize_to_layout!(fig)
    return fig
end

"Plot illustrative phase plot for overview figure"
function plot_illustrative(
    ;
    N=1.0, ρ=1.0, α=[0.5, 2.0], γ=1.0
    )
    #~ field copied from Symbolics output
    field(η,θ,α) = (
        ((N^3)*exp(θ)*ρ - (N^2)*exp(θ)*η*ρ - (N^2)*exp(-θ)*α*η - (N^2)*exp(-θ)*γ*η
         + 2.0N*exp(-θ)*α*(η^2) + N*exp(-θ)*γ*(η^2) - exp(-θ)*α*(η^3)) / (N^2),
        (-(N^2)*ρ + (N^2)*exp(θ)*ρ - (N^2)*exp(-θ)*α - (N^2)*exp(-θ)*γ +
         4.0N*exp(-θ)*α*η + 2.0N*exp(-θ)*γ*η - 3.0exp(-θ)*α*(η^2) +
         (N^2)*exp(θ)*exp(-θ)*α + (N^2)*exp(θ)*exp(-θ)*γ - 4.0N*exp(θ)*exp(-θ)*α*η -
         2.0N*exp(θ)*exp(-θ)*γ*η + 3.0exp(θ)*exp(-θ)*α*(η^2)) / (N^2)
    )
    #/ Create figure
    fig = Figure(; size=(128,278), backgroundcolor=:transparent, figure_padding=(0,2,2,0))

    #~ order param. plot
    axφ = Axis(
        fig[2,1], width=64, height=64, limits=(0.5,1.5,0,.8),
        yticks = [0.0,0.8], yminorticks=[0.4], yminorticksvisible=true,
        xticksize=2, yticksize=2, yminorticksize=1.5,
        xlabel=L"\textrm{control\;param.}\;\gamma", xlabelsize=10,
        ylabel=L"$\textrm{order\;param.}$\n$\varphi=1-\eta$", ylabelsize=10,
        xticklabelsvisible=false, yticklabelsvisible=false,
        xlabelpadding=0, ylabelpadding=0,
        xgridvisible=false, ygridvisible=false
    )
    #~ phase portraits 
    ax2 = Axis(
        fig[3,1], width=64, height=64, limits=(0,2,-0.5,1),        
        xlabel=L"\eta", ylabel=L"\theta", xlabelsize=10, ylabelsize=10,
        xlabelpadding=0, ylabelpadding=0,
        xticklabelsvisible=false, yticklabelsvisible=false,
        xticksvisible=false, yticksvisible=false
    )
    ax1 = Axis(
        fig[1,1], width=64, height=64, limits=(0,2,-0.5,1),
        xlabel=L"\eta", ylabel=L"\theta", xlabelsize=10, ylabelsize=10,
        xlabelpadding=0, ylabelpadding=0,
        # xlabelvisible=false, ylabelvisible=false,
        xticklabelsvisible=false, yticklabelsvisible=false,
        xticksvisible=false, yticksvisible=false
    )
    # ax1 = Axis(gl[1,1])  #~ 1st order transition plot
    # ax2 = Axis(gl[3,1])  #~ 2nd order transition plot
    textlabels = [L"2\textrm{nd\;order}", L"1\textrm{st\;order}"]
    axes = [ax2, ax1]
    for i in eachindex(α)
        f(η,θ) = Point2f(field(η,θ,α[i]))
        sp = streamplot!(
            axes[i], f,
            0.0..2.0, -0.5..1.0, color= x -> :lightgray, alpha=.7,
            arrow_size=4., linewidth=.25,
            density=.5, maxsteps=500, stepsize=0.001,
            gridsize=(32,32)
        )
        band!(axes[i], [1,2], [-0.5,-0.5], [1.0,1.0], color=(:red, 0.2))

        # / Add trivial zero-energy lines
        vlines!(axes[i], [1.0], color=:black, linewidth=1.2) 
        hlines!(axes[i], [0.0], color=:black, linewidth=1.2)
        
        #/ Add non-trivial zero-energy line
        ηmax = min(2.0, N*(α[i]+γ)/α[i] - 1e-7)
        ηplot = 0.0:0.01:ηmax
        θf(η) = @. log(η * (α[i]*(N-η) + γ*N)/(ρ*N^2))
        lp = lines!(
            axes[i], ηplot, θf(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )

        #/ Add arrows
        #~ field on trivial zero-energy lines with θ=0
        ηt = α[i] < ρ ? [0.5, 1.5] : [0.25, 0.75, 1.5]
        plot_arrowfield(axes[i], (x,y)->field(x,y,α[i]), ηt, zeros(length(ηt)), :black)
        θt = γ < ρ ? [-0.2, 0.5] : ((γ == ρ) ? [-0.3, 0.5] : [-0.3, 0.35])
        plot_arrowfield(axes[i], (x,y)->field(x,y,α[i]), ones(length(θt)), θt, :black)
        #~ field on non-trivial zero-energy line(s)
        ηnt = [0.3, 0.75, (1.0+ηmax)/2]
        θnt = θf.(ηnt)
        plot_arrowfield(axes[i], (x,y)->field(x,y,α[i]), ηnt, θnt, :rebeccapurple)
        
        #/ Scatter fixed points        
        if i == 1
            #~ fixed point F below critical transition
            F = [[1.0, 1.0], [0.0, θf(1.0)]]
            markers = [:circle, :star5]
            markersize = [6, 9]
            colors = [:mediumpurple1, :white]
            strokecolors = [:rebeccapurple, :firebrick2]            
            s = scatter!(
                axes[i], F[begin], F[end], marker=markers, markersize=markersize,
                color=colors, strokewidth=.9, strokecolor=strokecolors           
            )
        else i == 2
            # #~ fixed points at critical transition
            ηmin = ((α[i]+γ)*N - sqrt((α[i]+γ)^2*N^2 - 4*α[i]*N^3*ρ)) / (2*α[i]*N)
            F = [[1.0, ηmin], [0.0, 0.0]]
            markers = [:circle, :star5]
            markersize = [6,9]
            colors = [:mediumpurple1, :white]
            strokecolors = [:rebeccapurple, :firebrick2]
            labels = [L"\textrm{global fixed point}", L"\textrm{fixed point}"]
            s = scatter!(
                axes[i], F[begin], F[end], marker=markers, markersize=markersize,
                color=colors, strokewidth=.9, strokecolor=strokecolors
            )
            s2 = scatter!(axes[i], [-2.0,-2.0], marker=markers[begin], color=colors[begin],
                strokewidth=.9, strokecolor=strokecolors[begin]
            )
            s1 = scatter!(axes[i], [-2.0,-2.0], marker=markers[end], color=colors[end],
                strokewidth=.9, strokecolor=strokecolors[end]
            )
            legend = Legend(
                fig[0,1], [s1, s2], labels, labelsize=8, nbanks=2,
                framevisible=false, padding=3, rowgap=0, patchsize=(3,3)
            )
        end
        # else        
        #     # #~ fixed points above critical transition
        #     ηmin = ((α[i]+γ)*N - sqrt((α[i]+γ)^2*N^2 - 4*α[i]*N^3*ρ)) / (2*α[i]*N)
        #     ηplus = ((α[i]+γ)*N + sqrt((α[i]+γ)^2*N^2 - 4*α[i]*N^3*ρ)) / (2*α[i]*N)
        #     F = [[1.0, ηmin, ηplus, 1.0], [0.0, 0.0, 0.0, θf(1.0)]]
        #     markers = [:circle, :star5, :circle, :circle]
        #     markersize = [6,9,6,6]
        #     colors = [:mediumpurple1, :white, :mediumpurple1, :mediumpurple1]
        #     strokecolors = [:rebeccapurple, :firebrick2, :rebeccapurple, :rebeccapurple]
        # end
        s = scatter!(
            axes[i], F[begin], F[end], marker=markers, markersize=markersize,
            color=colors, strokewidth=.9, strokecolor=strokecolors           
        )
        #~ some text
        text!(
            axes[i], 0.88, 0.7, text=L"\textbf{\textrm{non}}", space=:relative, fontsize=8,
            align=(:center,:bottom), color=:firebrick2, rotation=π/2
        )
        text!(
            axes[i], 0.98, 0.7, text=L"\textbf{\textrm{physical}}", space=:relative,
            fontsize=8, align=(:center,:bottom), color=:firebrick2, rotation=π/2
        )
        text!(
            axes[i], 0.02, 0.98, text=textlabels[i], space=:relative,
            fontsize=7, align=(:left,:top)
        )
    end
    for i in reverse(eachindex(α))
    #/ Plot order param.
        #/ Define function for the order parameter (from theory, see below)
        ϕ(γ) = 1 - ((α[i]+γ)*N - sqrt((α[i]+γ)^2*N^2 - 4*α[i]*N^3*ρ)) / (2*α[i]*N)
        #~ specify γ values
        γv = ρ:0.01:2ρ
        αlabel = α[i] < ρ ? L"\alpha < \rho" : L"\alpha > \rho"
        color = α[i] < ρ ? :firebrick2 : :firebrick2
        linestyle = α[i] < ρ ? :solid : (:dash, :dense)

        #/ Plot lines
        lines!(
            axφ, vcat(0.0, ρ), vcat(0.0, 0.0), color=color, linewidth=1.2, label=αlabel,
            linestyle=linestyle
        )
        lines!(axφ, γv, ϕ.(γv), color=color, linewidth=1.2, linestyle=linestyle)
        #~ if 1st order transition, identify it more clearly    
        if α[i] > ρ
            vlines!(
                axφ, [1.0], ymin=0.0, ymax=0.6, color=:gray,
                linewidth=1., linestyle=(:dot,:dense)
            )
        end
        x = (α[i] > ρ) ? 1.1 : 1.075
        y = (α[i] > ρ) ? 0.205 : 0.57
        rotation = (α[i] > ρ) ? π/5 : π/12
        text!(
            axφ, x, y, text=textlabels[mod1(i+1,2)], align=(:left,:bottom), fontsize=7,
            color=:black, rotation=rotation
        )
    end
    
    # axislegend(
    #     axφ, position=:lt, labelsize=9, framevisible=false, rowgap=0,
    #     patchsize=(10,1), padding=0
    # )

    
    
    rowgap!(fig.layout, 5)
    rowgap!(fig.layout, 1, 2)
    resize_to_layout!(fig)
    return fig 
end

"Plot spatially extended predator-prey model (Lotka-Volterra)"
function plot_phase_sLV(; D=[0.3,0.8], c=1.0, k=0.4, ε=0.2, N=1.0)
    #/ Define figure
    width = .8*246
    fig = Figure(;
        size=(width,width/1.67), figure_padding=(1,2,2,1),
        backgroundcolor=:transparent
    )
    #/ Define axes
    axtitles = [L"D = D_1", L"D = D_2"]
    ax = [Axis(
        fig[1,i],
        title = axtitles[i], titlesize=10, titlegap=1,
        xlabel=L"\eta", xlabelpadding=0.0,
        ylabel=L"\theta", ylabelpadding=0.0,
        aspect=1,
        xlabelsize=12, ylabelsize=12, ylabelvisible=(i==1),
        xticks = [0.0, 1.0, 2.0], xticksize=2,
        xminorticksvisible=true, xminorticks=IntervalsBetween(2), xminorticksize=1,
        yticks=[-0.5, 0.0, 0.5], yticksize=2,
        xticklabelsize=7, yticklabelsize=7, yticklabelsvisible=(i==1),
        xgridvisible=false, ygridvisible=false,
        limits=(-.5,1.5,-.5,.5),
    ) for i in 1:2]

    #~ Specify critical points
    D1 = N - ε - ε/k
    D2 = N - ε/c
    @info "D" D1 D2

    for (i,d) in enumerate(D)
        #/ Add non-physical band(s)
        band!(ax[i], [1,2], [-0.5,-0.5], [1.0,1.0], color=(:red, 0.2))
        band!(ax[i], [-.5,0.], [-0.5,-0.5], [0.5,0.5], color=(:red, 0.2))
        #/ Add streamplot 
        #~ Specify the field (derived symbolically, see symbolics/LV.jl)
        field(η,θ) = d ≤ D1 ?
                    ((-(N^3)*(ε^3)*η + (N^3)*c*exp(θ)*(ε^2)*η - (N^3)*exp(-θ)*(ε^3)*η - (N^2)*(c^2)*exp(-θ)*ε*(η^2) - (N^2)*c*d*exp(θ)*(ε^2)*η - (N^2)*c*k*exp(-θ)*ε*(η^2) - (2//1)*(N^2)*c*exp(-θ)*(ε^2)*(η^2) + (N^2)*k*exp(-θ)*(ε^2)*(η^2) + N*(c^2)*d*exp(-θ)*ε*(η^2) + N*c*d*k*exp(-θ)*ε*(η^2) + N*c*k*exp(-2θ)*ε*(η^3) + (N^3)*exp(θ)*exp(-θ)*(ε^3)*η + (2//1)*(N^2)*(c^2)*exp(θ)*exp(-θ)*ε*(η^2) - (2//1)*(N^2)*c*k*(exp(-θ)^2)*ε*(η^2) - (2//1)*(N^2)*c*(exp(-θ)^2)*(ε^2)*(η^2) + (2//1)*(N^2)*k*(exp(-θ)^2)*(ε^2)*(η^2) - (2//1)*N*(c^2)*d*exp(θ)*exp(-θ)*ε*(η^2) + N*(c^2)*k*exp(-2θ)*exp(-θ)*(η^3) - N*(c^2)*k*(exp(-θ)^2)*(η^3) - N*(c^2)*(exp(-θ)^2)*ε*(η^3) + (2//1)*N*c*d*k*(exp(-θ)^2)*ε*(η^2) + (2//1)*N*c*k*exp(-2θ)*exp(-θ)*ε*(η^3) + N*c*k*(exp(-θ)^2)*ε*(η^3) - (c^2)*d*k*exp(-2θ)*exp(-θ)*(η^3) + (c^2)*d*k*(exp(-θ)^2)*(η^3) - (c^2)*k*(exp(-2θ)^2)*(η^4) + (c^2)*k*exp(-2θ)*exp(-θ)*(η^4) + (2//1)*(N^2)*c*k*exp(θ)*(exp(-θ)^2)*ε*(η^2) + (2//1)*(N^2)*c*exp(θ)*(exp(-θ)^2)*(ε^2)*(η^2) - (2//1)*(N^2)*k*exp(θ)*(exp(-θ)^2)*(ε^2)*(η^2) - N*(c^2)*k*exp(-2θ)*exp(θ)*exp(-θ)*(η^3) - (2//1)*N*(c^2)*k*(exp(-θ)^3)*(η^3) - N*(c^2)*(exp(-θ)^3)*ε*(η^3) - (2//1)*N*c*d*k*exp(θ)*(exp(-θ)^2)*ε*(η^2) - (2//1)*N*c*k*exp(-2θ)*exp(θ)*exp(-θ)*ε*(η^3) + (2//1)*N*c*k*(exp(-θ)^3)*ε*(η^3) + (c^2)*d*k*exp(-2θ)*exp(θ)*exp(-θ)*(η^3) + (2//1)*(c^2)*d*k*(exp(-θ)^3)*(η^3) + (c^2)*k*(exp(-2θ)^2)*exp(θ)*(η^4) + (3//1)*(c^2)*k*exp(-2θ)*(exp(-θ)^2)*(η^4) + (2//1)*N*(c^2)*k*exp(θ)*(exp(-θ)^3)*(η^3) + N*(c^2)*exp(θ)*(exp(-θ)^3)*ε*(η^3) - (2//1)*N*c*k*exp(θ)*(exp(-θ)^3)*ε*(η^3) - (2//1)*(c^2)*d*k*exp(θ)*(exp(-θ)^3)*(η^3) - (3//1)*(c^2)*k*exp(-2θ)*exp(θ)*(exp(-θ)^2)*(η^4)) / (N*((-N*ε - c*exp(-θ)*η)^2)), ((N^3)*c*(ε^2) - (N^2)*c*d*(ε^2) - (N^3)*c*exp(θ)*(ε^2) - (N^3)*exp(-θ)*(ε^3) + (N^2)*c*d*exp(θ)*(ε^2) + (N^3)*exp(θ)*exp(-θ)*(ε^3) - (2//1)*(N^2)*c*k*(exp(-θ)^2)*ε*η - (2//1)*(N^2)*c*(exp(-θ)^2)*(ε^2)*η + (2//1)*(N^2)*k*(exp(-θ)^2)*(ε^2)*η + N*(c^2)*k*exp(-2θ)*exp(-θ)*(η^2) + (2//1)*N*c*d*k*(exp(-θ)^2)*ε*η + (2//1)*N*c*k*exp(-2θ)*exp(-θ)*ε*(η^2) - (c^2)*d*k*exp(-2θ)*exp(-θ)*(η^2) - (c^2)*k*(exp(-2θ)^2)*(η^3) + (2//1)*(N^2)*c*k*exp(θ)*(exp(-θ)^2)*ε*η + (2//1)*(N^2)*c*exp(θ)*(exp(-θ)^2)*(ε^2)*η - (2//1)*(N^2)*k*exp(θ)*(exp(-θ)^2)*(ε^2)*η - N*(c^2)*k*exp(-2θ)*exp(θ)*exp(-θ)*(η^2) - (2//1)*N*(c^2)*k*(exp(-θ)^3)*(η^2) - N*(c^2)*(exp(-θ)^3)*ε*(η^2) - (2//1)*N*c*d*k*exp(θ)*(exp(-θ)^2)*ε*η - (2//1)*N*c*k*exp(-2θ)*exp(θ)*exp(-θ)*ε*(η^2) + (2//1)*N*c*k*(exp(-θ)^3)*ε*(η^2) + (c^2)*d*k*exp(-2θ)*exp(θ)*exp(-θ)*(η^2) + (2//1)*(c^2)*d*k*(exp(-θ)^3)*(η^2) + (c^2)*k*(exp(-2θ)^2)*exp(θ)*(η^3) + (3//1)*(c^2)*k*exp(-2θ)*(exp(-θ)^2)*(η^3) + (2//1)*N*(c^2)*k*exp(θ)*(exp(-θ)^3)*(η^2) + N*(c^2)*exp(θ)*(exp(-θ)^3)*ε*(η^2) - (2//1)*N*c*k*exp(θ)*(exp(-θ)^3)*ε*(η^2) - (2//1)*(c^2)*d*k*exp(θ)*(exp(-θ)^3)*(η^2) - (3//1)*(c^2)*k*exp(-2θ)*exp(θ)*(exp(-θ)^2)*(η^3)) / (N*((-N*ε - c*exp(-θ)*η)^2))) :
                     ((-N*ε*η + N*c*exp(θ)*η - N*exp(-θ)*ε*η - c*d*exp(θ)*η - c*exp(-θ)*(η^2) + N*exp(θ)*exp(-θ)*ε*η) / N, (N*c - c*d - N*c*exp(θ) - N*exp(-θ)*ε + c*d*exp(θ) - 2c*exp(-θ)*η + N*exp(θ)*exp(-θ)*ε + 2c*exp(θ)*exp(-θ)*η) / N)
                     
                     # (
                     #     (c*exp(θ)*ε*η - c*exp(-θ)*ε*η - k*exp(θ)*ε*η -
                     #      d*c*k*exp(θ)*η + N*c*k*exp(θ)*η - c*k*exp(θ)*(η^2)) / c,
                     #     (c*ε - k*ε - d*c*k + N*c*k - 2c*k*η - c*exp(θ)*ε - c*exp(-θ)*ε
                     #      + k*exp(θ)*ε + d*c*k*exp(θ) - N*c*k*exp(θ) + 2c*k*exp(θ)*η +
                     #      c*exp(θ)*exp(-θ)*ε) / c) :
                     # (
                     #     -((d - N + η)*c + exp(-θ)*ε)*exp(θ)*η + (1 - exp(θ))*exp(-θ)*ε*η,
                     #     ((d - N + η)*c + exp(-θ)*ε)*(1 - exp(θ)) + c*(1 - exp(θ))*η
                     # )
        f(η,θ) = Point2f(field(η,θ))
        sp = streamplot!(
            ax[i], f,
            -0.5..1.5, -0.5..0.5, color= x -> :gray, alpha=.8,
            arrow_size=4., linewidth=.3,
            density=.5, maxsteps=512, stepsize=0.001,
            gridsize=(32,32)
        )

        #/ Add trivial zero-energy lines
        vlines!(ax[i], [0.0], color=:black, linewidth=1.0) 
        hlines!(ax[i], [0.0], color=:black, linewidth=1.0)
        #/ Add non-trivial zero-energy line
        θf(η) = d ≤ D1 ?
                log(-ε^2 + sqrt(ε^4 - 4*(c*d*ε-c*ε)*(c*k*η-c*d*k*η+c*ε*η-k*ε*η-c*k*η^2)) / (2*(-c*ε + c*d*ε))) :
               # -log(N - k/c + k/ε*(N-d-η)) :
               -log(c/ε * (N-d-η))
        θf_incorrect(η) = d ≤ D1 ?
               -log(c/ε * (N-d-η)) :
               -log(N - k/c + k/ε*(N-d-η))
        ηmax = d ≤ D1 ? N - d : 0.1 #N - d - ε/c + N*ε/k
        # ηmax_incorrect = d < D1 ? 0.5 : 0.5
        # ηmin_incorrect = d < D1 ? 0.4 : -0.375
        # ηplot_incorrect = ηmin_incorrect:0.01:ηmax_incorrect
        # lines!(
        #     ax[i], ηplot_incorrect, θf_incorrect.(ηplot_incorrect),
        #     color=:black, linewidth=.8, linestyle=(:dot,:dense)
        # )
        ηplot = -0.2:0.01:ηmax
        @info "hmm" θf.(ηplot)
        lines!(
            ax[i], ηplot, θf.(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )
        # ηplot = -0.5:0.01:0.5

        F = [[0.0, N - d - ε/c], [0.0,0.0]]
        markers = [:circle, :star5]
        markersize = [6,9]
        colors = [:mediumpurple1, :white]
        strokecolors = [:rebeccapurple, :firebrick2]
        s = scatter!(
            ax[i], F[begin], F[end], marker=markers, markersize=markersize,
            color=colors, strokewidth=.9, strokecolor=strokecolors
        )
        #/ Add arrows
        #~ field on trivial zero-energy lines with θ=0
        ηt = d ≤ D1 ? [-0.25,0.25,1.0] : [-0.25, 0.5]
        plot_arrowfield(ax[i], (x,y)->field(x,y), ηt, zeros(length(ηt)), :black)
        θt = [-0.25, 0.25]
        plot_arrowfield(ax[i], (x,y)->field(x,y), 0.0.*ones(length(θt)), θt, :black)
        ηnt = d ≤ D1 ? [0.3, 0.65] : [-0.1, 0.055]
        θnt = θf.(ηnt)
        plot_arrowfield(ax[i], (x,y)->field(x,y), ηnt, θnt, :rebeccapurple)
    end

    #~ Reduce some white-space
    colgap!(fig.layout, 6)
    resize_to_layout!(fig)
    return fig 
end

"Plot phase-plot for tax-evasion model"
function plot_phase_taxevasion(; N=1.0, α=1.0, β=0.5, δ=0.5)
    #/ Define figure
    width = .8*246
    fig = Figure(;
        size=(width,width/1.67), figure_padding=(1,2,2,1),
        backgroundcolor=:transparent
    )
    #/ Define axes
    axtitles = [L"\epsilon < \epsilon_c", L"\epsilon > \epsilon_c"]
    ax = [Axis(
        fig[1,i],
        title = axtitles[i], titlesize=10, titlegap=1,
        xlabel=L"\eta", xlabelpadding=0.0,
        ylabel=L"\theta", ylabelpadding=0.0,
        aspect=1,
        xlabelsize=12, ylabelsize=12, ylabelvisible=(i==1),
        xticks = [0.0, 1.0, 2.0], xticksize=2,
        xminorticksvisible=true, xminorticks=IntervalsBetween(2), xminorticksize=1,
        yticks=[-0.5, 0.0, 0.5], yticksize=2,
        xticklabelsize=7, yticklabelsize=7, yticklabelsvisible=(i==1),
        xgridvisible=false, ygridvisible=false,
        limits=(0.0,1.5,-1.,1.),
    ) for i in 1:2]

    #/ Compute εc
    εc = α*β/δ
    εv = [εc/3, 3*εc]

    for (i,ε) in enumerate(εv)    
        #/ Compute λc, as we want to plot at λ=λc
        λ = (α + ε)*(β + δ) / α
    
        #/ Add non-physical band(s)
        band!(ax[i], [1,2], [-1.0,-1.0], [1.0,1.0], color=(:red, 0.2))
        #/ Add streamplot 
        #~ Specify the field (derived symbolically, see symbolics/tax-evasion.jl)
        field(η,θ) = (
            ((-exp(θ)*(α^2)*η*λ - exp(θ)*α*β*δ*η - exp(θ)*α*β*η*λ + (exp(θ)^3)*(α^2)*β + (exp(θ)^3)*α*(β^2) + (2//1)*(exp(θ)^2)*α*β*δ*η + exp(θ)*(α^2)*(η^2)*λ + exp(θ)*α*β*δ*(η^2) + exp(θ)*α*β*(η^2)*λ - exp(θ)*α*(δ^2)*(η^2) + exp(θ)*α*δ*ε*(η^2) - exp(θ)*α*δ*(η^2)*λ - (exp(θ)^3)*(α^2)*β*η + (exp(θ)^3)*(α^2)*δ*η - (exp(θ)^3)*α*(β^2)*η + (exp(θ)^3)*α*β*δ*η + (exp(θ)^3)*α*β*ε*η + (exp(θ)^3)*(β^2)*ε*η - (2//1)*(exp(θ)^2)*α*β*δ*(η^2) + (2//1)*(exp(θ)^2)*α*(δ^2)*(η^2) + (2//1)*(exp(θ)^2)*β*δ*ε*(η^2) + exp(θ)*α*(δ^2)*(η^3) - exp(θ)*α*δ*ε*(η^3) + exp(θ)*α*δ*(η^3)*λ + exp(θ)*(δ^2)*ε*(η^3) - (exp(θ)^3)*(α^2)*δ*(η^2) - (exp(θ)^3)*α*β*δ*(η^2) - (exp(θ)^3)*α*β*ε*(η^2) - (exp(θ)^3)*(β^2)*ε*(η^2) - (2//1)*(exp(θ)^2)*α*(δ^2)*(η^3) - (2//1)*(exp(θ)^2)*β*δ*ε*(η^3) - exp(θ)*(δ^2)*ε*(η^4)) / ((exp(θ)*α + exp(θ)*β + δ*η)^2),
             (-exp(θ)*(α^2)*λ - exp(θ)*α*β*δ - exp(θ)*α*β*λ - (exp(θ)^2)*(α^2)*β + (exp(θ)^2)*(α^2)*δ + (exp(θ)^2)*(α^2)*λ - (exp(θ)^2)*α*(β^2) + (2//1)*(exp(θ)^2)*α*β*δ + (exp(θ)^2)*α*β*ε + (exp(θ)^2)*α*β*λ + (exp(θ)^2)*(β^2)*ε + (2//1)*exp(θ)*(α^2)*η*λ + (2//1)*exp(θ)*α*β*η*λ + (2//1)*exp(θ)*α*δ*ε*η + (2//1)*exp(θ)*β*δ*ε*η + α*δ*(η^2)*λ + (δ^2)*ε*(η^2) + (exp(θ)^3)*(α^2)*β - (exp(θ)^3)*(α^2)*δ + (exp(θ)^3)*α*(β^2) - (exp(θ)^3)*α*β*δ - (exp(θ)^3)*α*β*ε - (exp(θ)^3)*(β^2)*ε - (2//1)*(exp(θ)^2)*(α^2)*δ*η - (2//1)*(exp(θ)^2)*(α^2)*η*λ - (2//1)*(exp(θ)^2)*α*β*δ*η - (2//1)*(exp(θ)^2)*α*β*ε*η - (2//1)*(exp(θ)^2)*α*β*η*λ - (2//1)*(exp(θ)^2)*α*δ*ε*η - (2//1)*(exp(θ)^2)*(β^2)*ε*η - (2//1)*(exp(θ)^2)*β*δ*ε*η - exp(θ)*α*(δ^2)*(η^2) - (3//1)*exp(θ)*α*δ*ε*(η^2) - exp(θ)*α*δ*(η^2)*λ - (4//1)*exp(θ)*β*δ*ε*(η^2) - exp(θ)*(δ^2)*ε*(η^2) - (2//1)*(δ^2)*ε*(η^3) + (2//1)*(exp(θ)^3)*(α^2)*δ*η + (2//1)*(exp(θ)^3)*α*β*δ*η + (2//1)*(exp(θ)^3)*α*β*ε*η + (2//1)*(exp(θ)^3)*(β^2)*ε*η + (exp(θ)^2)*α*(δ^2)*(η^2) + (3//1)*(exp(θ)^2)*α*δ*ε*(η^2) + (4//1)*(exp(θ)^2)*β*δ*ε*(η^2) + (2//1)*exp(θ)*(δ^2)*ε*(η^3)) / ((exp(θ)*α + exp(θ)*β + δ*η)^2))
        )
        f(η,θ) = Point2f(field(η,θ))
        sp = streamplot!(
            ax[i], f,
            0.0..1.5, -1.0..1.0, color= x -> :gray, alpha=.8,
            arrow_size=4., linewidth=.3,
            density=.5, maxsteps=512, stepsize=0.001,
            gridsize=(32,32)
        )

        #/ Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:black, linewidth=1.0) 
        hlines!(ax[i], [0.0], color=:black, linewidth=1.0)
        #/ Add non-trivial zero-energy line
        θf(η) = log( (α*λ*η - δ*ε*η^2) / (α*β + α*δ*η + β*ε*η) )
        ηmin = 0.0
        ηmax = 1.5
        
        ηplot = ηmin:0.01:ηmax
        lines!(
            ax[i], ηplot, θf.(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )
        
        _η = (α*(λ-δ)-β*ε-δ*ε*sqrt(β^2*ε^2 + α^2*(λ-δ)^2 - 2*α*β*ε*(λ+δ))/(δ*ε)) / (2*δ*ε)
        F = ε<εc ? [[-1.0, 1.0], [-1.0, 0.0]] : [[1.0,_η], [0.0,0.0]]
        markers = [:circle, :star5]
        markersize = [6,9]
        colors = [:mediumpurple1, :white]
        strokecolors = [:rebeccapurple, :firebrick2]
        s = scatter!(
            ax[i], F[begin], F[end], marker=markers, markersize=markersize,
            color=colors, strokewidth=.9, strokecolor=strokecolors
        )
        #/ Add arrows
        #~ field on trivial zero-energy lines with θ=0 or η=1
        ηt = ε < εc ? [0.5,1.25] : [-0.25, 0.6]
        plot_arrowfield(ax[i], (x,y)->field(x,y), ηt, zeros(length(ηt)), :black)
        θt = [-0.5, 0.5]
        plot_arrowfield(ax[i], (x,y)->field(x,y), ones(length(θt)), θt, :black)
        #~ field on non-trivial zero-energy lines
        ηnt = ε < εc ? [0.5, 1.35] : [0.1, 0.675, 1.25]
        θnt = θf.(ηnt)
        plot_arrowfield(ax[i], (x,y)->field(x,y), ηnt, θnt, :rebeccapurple)
    end

    #~ Reduce some white-space
    colgap!(fig.layout, 6)
    resize_to_layout!(fig)
    return fig 
end

"Plot η-value (occupation no.) of global fixed point F vs γ, for both α"
function plot_occupation_taxevasion(; N=1.0, α=1.0, β=0.5, δ=0.5)
    #/ Create figures
    width = .8*246
    fig = Figure(
        ; size=(width,.6*width), figure_padding=(1,6,1,3),
        backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1],
        xlabel=L"\lambda",
        ylabel=L"\varphi",
        xlabelsize=11, ylabelsize=11,
        xlabelpadding=1.0, ylabelpadding=3.0,
        yticks=[0.0, 1.0], xticks=[0.0, 2.0, 4.0],
        xticklabelsize=8, yticklabelsize=8,
        xminorticksvisible=true, xminorticks=IntervalsBetween(4),
        yminorticksvisible=true, yminorticks=IntervalsBetween(4),
        limits=(0.0, 6.0, 0., 1.01),
        xgridvisible=false, ygridvisible=false
    )
    #/ Define function for the order parameter (from theory, see below)
    ϕ(λ,ε) = 1 - (α*(λ-δ) - β*ε - δ*ε*sqrt(β^2*ε^2 + α^2*(λ-δ)^2 - 2*α*β*ε*(λ+δ))/(δ*ε)) / (2*δ*ε)
    #~ specify γ values
    # αlabel = α < ρ ? L"\alpha < \rho" : L"\alpha > \rho"
    # color = α < ρ ? :firebrick2 : :firebrick2
    colors = [:black, :firebrick2]
    εlabels = [L"\epsilon\;<\;\epsilon_c", L"\epsilon\;>\;\epsilon_c"]

    #/ Plot lines
    # lines!(ax, vcat(0.0, ρ), vcat(0.0, 0.0), color=color, linewidth=1.2, label=αlabel)
    εc = α*β/δ
    εv = [εc/2, 3*εc]
    for i in eachindex(εv)
        λc = (α + εv[i])*(β + δ) / α
        @info λc
        λv = λc:0.01:6.0
    
        lines!(ax, λv, ϕ.(λv,εv[i]), color=colors[i], linewidth=1.2, label=εlabels[i])
        #~ if 1st order transition, identify it more clearly
        if εv[i] > εc
            vlines!(
                ax, [λc], ymin=0.0, ymax = ϕ(λv[begin], εv[i]), color=colors[i],
                linewidth=1., linestyle=:dot
            )
        end
    end

    axislegend(
        ax, position=:lt, labelsize=11, framevisible=false, rowgap=0,
        patchsize=(5,1), padding=0
    )
    return fig 
end

end # module PhasePlotter
#/ End module


