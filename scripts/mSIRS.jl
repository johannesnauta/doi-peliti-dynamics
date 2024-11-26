#= Module for simulating dynamics of modified SIRS model
 * includes functions that numerically integrate the reduced dynamics
=#
#/ Start module
module SIRSReDynamics

#/ Packages
using Catalyst
using DifferentialEquations
using FHist
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

#################
### FUNCTIONS ###
"Catalyst ReactionSystem"
function get_ReactionSystem()
    rn = @reaction_network begin
        γ, S + I --> 2I
        α, S + 2I --> 3I
        ρ, I --> R
        σ, R --> S
    end
    return rn    
end

function get_SDEProblem(
    rs::ReactionSystem;
    γ=1.0, α=0.4, ρ=0.4, σ=0.1, N=1e3, tspan=(0.0,1e4)
    ) 
    u = [0.95*N, 0.05*N, 0.0]
    umap = (:S => u[1], :I => u[2], :R => u[3])
    pmap = (:γ => γ/N, :α => α/N^2, :ρ => ρ, :σ => σ)
    sprob = SDEProblem(rs, umap, tspan, pmap)
    return sprob
end

function solve_SDE(sprob::SDEProblem; dt=1e-3)
    return solve(sprob, EM(), dt=dt)
end

function get_hist(sol::RODESolution; idx=1, tstart=2e2, bins=2e2:5e0:6e2)
    #~ Collect data
    x = reduce(hcat, sol.u)[idx,findfirst(sol.t .> tstart):end]
    fh = FHist.Hist1D(x; binedges=bins)
    return fh |> normalize
end

"Get reduced SDEProblem"
function get_rSDEProblem(; γ=1.0, α=0.4, ρ=0.4, σ=0.1, N=1e3, tspan=(0.0,1e4))
    #~ Define reduced SDEProblem
    @variables η(t)
    eq = [D(η) ~ (σ*(N-η)*(ρ*N^2*(ρ+σ) - N*η*(α*σ + γ*(ρ+σ)) + α*σ*η^2)) / (N^2*(ρ+σ)^2)]
    neq = [sqrt((σ*(N-η)*(ρ*N^2*(ρ+σ) + N*η*(α*σ + γ*(ρ+σ)) - α*σ*η^2)) / (N^2*(ρ+σ)^2))]
    @named ode = ODESystem(eq, t)
    @named sde = SDESystem(ode, neq)
    sprob = SDEProblem(complete(sde), [η => 0.95*N], tspan)
    return sprob
end

"ODESystem of mean-field equations"
function get_ODESystem()
    #/ Specify parameters and variables
    @parameters γ α ρ σ
    @variables S(t) I(t) R(t)
    eqns = [
        D(S) ~ σ*R - γ/(S+I+R)*S*I - α/(S+I+R)^2*S*I^2,
        D(I) ~ γ/(S+I+R)*S*I + α/(S+I+R)^2*S*I^2 - ρ*I,
        D(R) ~ ρ*I - σ*R
    ]
    @named sys = ODESystem(eqns, t, [S,I,R], [γ,α,ρ,σ])
    return complete(sys)
end

"SDESystem of ODESystem, with noise strength ξv"
function get_SDESystem(sys::ODESystem; ξv = 0.1)
    noiseqns = [ξv*sqrt(sys.S*sys.I), -ξv*sqrt(sys.S*sys.I), 0]
    @named sde = SDESystem(sys, noiseqns)
    return complete(sde)
end

"SDEProblem from SDESystem"
function get_SDEProblem(sys::SDESystem; γ=1.0, α=0.0, ρ=0.5, σ=0.1, N=1e3, tspan=(0.0,1e4))
    pmap = (sys.γ => γ, sys.α => α, sys.ρ => ρ, sys.σ => σ)
    u = [N / 2, N / 3]
    append!(u, N - sum(u))
    umap = (sys.S => u[1], sys.I => u[2], sys.R => u[3])
    prob = SDEProblem(sys, umap, tspan, pmap)
    return prob
end


end # module SIRSReDynamics
#/ End module
