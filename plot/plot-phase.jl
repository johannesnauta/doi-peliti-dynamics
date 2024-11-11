#/ Start module
module PhasePlotter

using CairoMakie
using LaTeXStrings

function plot_phase_SIR(;
    N::Float64=1.0, ρ::Float64=1.0, α::Float64=1.4, γ=[2.0/3.0,1.0,2.0]
    )
    f(η,θ,γ) = Point2f(
        (N-η) * (α*η*exp(-θ)/N + γ*η*exp(-θ)/N - ρ*exp(θ)),
        (N-η) * (1 - exp(θ)) * ( (α+γ)/N*exp(-θ) - 2*α*η/N^2*exp(-θ) )
        - (1 - exp(θ)) * (α*η/N^2*exp(-θ)*(N-η) + γ*η/N*exp(-θ) - ρ)
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
        limits=(0,2,-1,1)
    ) for i in 1:length(γ)]

    
    for i in eachindex(γ)
        #/ Add streamplot
        _f(η,θ) = f(η,θ,γ[i])
        sp = streamplot!(
            ax[i], _f,
            0.0..2.0, -1.0..1.0, color= x -> :gray,
            arrow_size=4., linewidth=.35,
            density=.5, maxsteps=500, stepsize=0.001,
            gridsize=(32,32)
        )
        # / Add trivial zero-energy lines
        vlines!(ax[i], [1.0], color=:rebeccapurple, linestyle=(:dash,:dense), linewidth=1.2) 
        hlines!(ax[i], [0.0], color=:darkorange, linewidth=1.2)
        
        #/ Add non-trivial zero-energy line
        ηplot = 0.01:0.01:1.5
        θf(η) = @. -log(ρ / (α*η/N^2 * (N-η) + γ[i]*η/N))
        lp = lines!(
            ax[i], ηplot, θf(ηplot),
            color=:rebeccapurple, linewidth=1.2, linestyle=(:dash,:dense)
        )
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
            0.0..2.0, -1.0..1.0, color= x -> :gray,
            arrow_size=4., linewidth=.35,
            density=.5, maxsteps=500, stepsize=0.001,
            gridsize=(32,32)
        )
        for θarrow in [-0.25, 0.25]
            spθ = streamplot!(
                ax[i], _f,
                [1.0], [θarrow], color=(x->:rebeccapurple), arrow_size=6., linewidth=.4
            )
        end
        for ηarrow in [0.75, 1.25]
            spη = streamplot!(
                ax[i], _f,
                [ηarrow], [0.0], color=(x->:darkorange), arrow_size=6., linewidth=.4
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
            spnt = streamplot!(
                ax[i], _f,
                [η_nontrivial[i]], [θf(η_nontrivial[i])],
                color=x->:rebeccapurple, arrow_size=6., linewidth=.4
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
    colgap!(fig.layout, 8) #rowgap!(fig.layout, 0)
    return fig
end

end # module PhasePlotter
#/ End module


