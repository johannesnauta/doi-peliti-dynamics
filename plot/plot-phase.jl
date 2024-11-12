#/ Start module
module PhasePlotter

using CairoMakie
using LaTeXStrings


"Plot some arrowheads at all (x,y) with the direction given by the field"
function plot_arrowfield(ax, field, x, y, color)
    derivatives = field.(x, y)
    dx, dy = first.(derivatives), last.(derivatives)
    arrows!(
        ax, x, y, dx, dy, color=color,
        linewidth=0., align=:center, arrowsize=7
    )
end

function plot_phase_SIR(;
    N::Float64=1.0, ρ::Float64=1.0, α::Float64=2.0, γ=[2.0/3.0,1.0,2.0]
    )    
    field(η,θ,γ) = (
        -exp(-θ)*(N-η) * (-N*(α+γ)*η + α*η^2 + exp(2*θ)*N^2*ρ) / N^2,
        -exp(-θ)*(-1+exp(θ)) * (-2*N*(2*α + γ)*η + 3*α*η^2 + N^2 * (α+γ + exp(θ)*ρ)) / N^2
    )
    f(η,θ,γ) = Point2f(field(η,θ,γ))
    
    width = .8*2*246
    fig = Figure(; size=(width,width/3), figure_padding=(2,8,2,2))
    axtitles = [L"\gamma < \rho", L"\gamma = \rho", L"\gamma > \rho"]
    
    ax = [Axis(
        fig[1,i], title=axtitles[i], titlesize=11,
        xlabel=L"\eta",
        ylabel=L"\theta",
        xlabelsize=12, ylabelsize=12,
        ylabelvisible=(i==1),            
        xticklabelsize=7, yticklabelsize=7,
        xgridvisible=false, ygridvisible=false,
        limits=(0,2,-.5,1)
    ) for i in 1:length(γ)]
    
    
    for i in eachindex(γ)
        band!(ax[i], [-1,2], [-0.5,-0.5], [0.0,0.0], color=(:red, 0.2))
        #/ Add streamplot
        _f(η,θ) = f(η,θ,γ[i])
        sp = streamplot!(
            ax[i], _f,
            0.0..2.0, -0.5..1.0,
            # colormap=:magma,
            color= x -> :gray, alpha=0.8,
            arrow_size=3.5, linewidth=.25,
            density=.3, maxsteps=256, stepsize=0.001,
            gridsize=(50,50)
        )
        # / Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:rebeccapurple, linestyle=(:dash,:dense), linewidth=1.2) 
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
        ηt = γ[i] < ρ ? [0.5, 1.5] : [0.75, 1.5]
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ηt, zeros(length(ηt)), :black)
        θt = γ[i] ≤ ρ ? [-0.5, 0.5] : [-0.5, 0.3]
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ones(length(θt)), θt, :rebeccapurple)
        #~ field on non-trivial zero-energy line(s)
        ηnt = γ[i] ≤ ρ ? [0.3, 0.7, (1.0+ηmax)/2] : [0.3, 0.6, (1.0+ηmax)/2]
        θnt = θf.(ηnt)        
        plot_arrowfield(ax[i], (x,y)->field(x,y,γ[i]), ηnt, θnt, :rebeccapurple)
        
        #/ Scatter fixed points
        if i == 1
            s = scatter!(
                ax[i], [1], [0], markersize=6, marker=:circle,
                color=:white, strokewidth=.9, strokecolor=:firebrick2
            )
        elseif i == 2
            scatter!(
                ax[i], [((α+γ[i])*N - sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)], [0.0],
                marker=:star5, color=:white, markersize=8.5,
                strokewidth=1., strokecolor=:firebrick2
            )
            s = scatter!(
                ax[i], [1], [0],
                markersize=6, marker=:circle, color=:white,
                strokewidth=.9, strokecolor=:rebeccapurple
            )
        else
            scatter!(
                ax[i], [((α+γ[i])*N - sqrt((α+γ[i])^2*N^2 - 4*α*N^3*ρ)) / (2*α*N)], [0.0],
                marker=:star5, color=:white, markersize=8.5,
                strokewidth=.8, strokecolor=:firebrick2
            )
            s = scatter!(
                ax[i], [1, 1], [0, θf(1.0)],
                markersize=6, marker=:circle, color=:white,
                strokewidth=.9, strokecolor=:rebeccapurple
            ) 
        end
        
    end
    return fig
end

function plot_phase_SIS(; N::Float64=1.0, ρ::Float64=1.0, γ = [2.0/3.0,1.0,2.0])
    f(η,θ,γ) = Point2f(
        (N-η) * (γ*η*exp(-θ)/N - ρ*exp(θ)),
        (1 - exp(θ)) * (γ*exp(-θ) - 2*γ*η/N*exp(-θ) + ρ)
    )
    width = .8*2*246
    fig = Figure(; size=(width,width/3), figure_padding=(2,8,2,2))
    axtitles = [L"\gamma < \rho", L"\gamma = \rho", L"\gamma > \rho"]
    
    ax = [Axis(
        fig[1,i], title=axtitles[i], titlesize=11,
        xlabel=L"\eta",
        ylabel=L"\theta",
        xlabelsize=12, ylabelsize=12,
        ylabelvisible=(i==1),            
        xticklabelsize=7, yticklabelsize=7,
        xgridvisible=false, ygridvisible=false,
        limits=(0,2,-1,1)
    ) for i in 1:length(γ)]

    strokecolors = [
        [:firebrick2,:black], [:firebrick2,:firebrick2], [:black,:firebrick2]
    ]
    
    for i in eachindex(γ)
        #/ Add streamplot
        _f(η,θ) = f(η,θ,γ[i])
        sp = streamplot!(
            ax[i], _f,
            0.0..2.0, -1.0..1.0, color= x -> :gray, alpha=.8,
            arrow_size=4., linewidth=.25,
            density=.5, maxsteps=500, stepsize=0.001,
            gridsize=(32,32)
        )
        for θarrow in [-0.25, 0.25]
            dη, dθ = _f(1.0, θarrow)
            arrθ = arrows!(
                ax[i], [1.0], [θarrow], [dη], [dθ], color=:rebeccapurple,
                linewidth=0., align=:center, arrowsize=7
            )
        end
        for ηarrow in [0.75, 1.25]
            dη, dθ = _f(ηarrow, 0.0)
            arrθ = arrows!(
                ax[i], [ηarrow], [0.0], [dη], [dθ], color=:darkorange,
                linewidth=0., align=:center, arrowsize=7
            )
        end
        # / Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:rebeccapurple, linestyle=(:dash,:dense), linewidth=1.2) 
        hlines!(ax[i], [0.0], color=:darkorange, linewidth=1.2)
        
        #/ Add non-trivial zero-energy line        
        ηplot = 0.01:0.01:2.0
        θf(η) = @. log(η * γ[i] / (ρ * N))
        lp = lines!(
            ax[i], ηplot, θf(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )
        
        #/ Scatter fixed points
        η_nontrivial = [1.25, nothing, 0.75]
        if i!=2
            dη, dθ = _f(η_nontrivial[i], θf(η_nontrivial[i]))
            arrθ = arrows!(
                ax[i], [η_nontrivial[i]], [θf(η_nontrivial[i])], [dη], [dθ],
                color=:rebeccapurple, linewidth=0., align=:center, arrowsize=7
            )
            snontrivial = scatter!(
                ax[i], [1.0], [log(γ[i]/ρ)], markersize=6,
                color=:white, strokewidth=.9, strokecolor=:rebeccapurple
            )
            s = scatter!(
                ax[i], [1,ρ*N/γ[i]], [0,0], markersize=6,
                color=:white, strokewidth=.9, strokecolor=[:black,:firebrick2]
            )
        else
            s = scatter!(
                ax[i], [1], [0], markersize=6,
                color=:white, strokewidth=.9, strokecolor=:firebrick2
            )
        end
        
    end
    colgap!(fig.layout, 8)
    return fig
end


end # module PhasePlotter
#/ End module


