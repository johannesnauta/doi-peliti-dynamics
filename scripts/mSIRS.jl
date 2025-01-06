#= Module for simulating dynamics of modified SIRS model
 * includes functions that numerically integrate the reduced dynamics
 * if not specified, closed-form derivations of the reduced dynamics have been
   performed using Mathematica
=#

#/ Start module
module SIRSReDynamics

#/ Packages
using Catalyst
using DifferentialEquations
using FHist
using JLD2
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

#################
### FUNCTIONS ###
function solve_SIR(;
    N=1e4, γ=1.0, α=0.2, ρ=0.4, σ=1.0, dt=1e-2, tspan=(0.0,1e5), meanfield=false
)
    rn = get_ReactionSystem()
    umap = (:S => 0.95*N, :I => 0.02*N, :R => 0.03*N)
    pmap = (:γ => γ/N, :α => α/N^2, :ρ => ρ, :σ => σ)
    sprob = SDEProblem(rn, umap, tspan, pmap)
    ssol = solve_SDE(sprob, dt=dt)
    #/ Solve mean-field eqns by solving an ODEProblem as well
    if meanfield
        oprob = ODEProblem(rn, umap, tspan, pmap)
        osol = solve_ODE(oprob)
        return ssol, osol
    end
    return ssol
end

function solve_rSIR(; N=1e4, γ=1.0, α=0.2, ρ=0.4, σ=1.0, tspan=(0.0,1e5))
    rsprob = get_rSDEProblem(; N=N, γ=γ, α=α, ρ=ρ, σ=σ, tspan=tspan)
    rssol = solve_SDE(rsprob, dt=1e-2)
    return rssol
end 

"Catalyst ReactionSystem"
function get_ReactionSystem()
    rn = @reaction_network begin
        γ, S + I --> 2I
        α*I, S + I --> 2I
        ρ, I --> R
        σ, R --> S
    end
    return rn    
end

function get_SDEProblem(
    rs::ReactionSystem;
    γ=1.0, α=0.0, ρ=0.1, σ=0.1, N=1e3, tspan=(0.0,5e2)
    ) 
    u = [0.95*N, 0.05*N, 0.0]
    umap = (:S => u[1], :I => u[2], :R => u[3])
    pmap = (:γ => γ/N, :α => α/N^2, :ρ => ρ, :σ => σ)
    sprob = SDEProblem(rs, umap, tspan, pmap)
    return sprob
end

"Get reduced SDEProblem"
function get_rSDEProblem(; γ=0.1, α=0.1, ρ=0.01, σ=0.01, N=1e3, tspan=(0.0,1e4))
    #~ Define reduced SDEProblem
    @variables η(t)
    eq = [
        D(η) ~ (N-η)*σ*(α*η^2*σ + N^2*ρ*(ρ+σ) - N*η*(α*σ + γ*(ρ+σ)))/(N^2*(ρ+σ)^2)
    ]
    neq = [
        sqrt((N-η)*σ*(-α*η^2*σ + N^2*ρ*(ρ+σ) + N*η*(α*σ + γ*(ρ+σ)))/(N^2*(ρ+σ)^2))
    ]    
    @named ode = ODESystem(eq, t)
    @named sde = SDESystem(ode, neq)
    sprob = SDEProblem(complete(sde), [η => 0.95*N], tspan)
    return sprob
end

"Solve trajectories for the full system (stochastic & mean-field) and reduced system"
function solve_trajectories(; tspan=(0.0, .5e2), dt=1e-2, ddir="../data/trajectories/")
    ssol, osol = solve_SIR(; tspan=tspan, dt=dt, meanfield=true)
    rsol = solve_rSIR(; tspan=tspan)
    #/ Extract trajectories for i=1 (S)
    xsde = reduce(hcat, ssol.u)[1,:]
    tsde = ssol.t
    xode = reduce(hcat, osol.u)[1,:]
    tode = osol.t
    xr   = reduce(hcat, rsol.u)[1,:]
    tr   = rsol.t
    #~ Save
    fname = ddir*"mSIRStrajectories.jld2"
    jldsave(fname; xsde=xsde, tsde=tsde, xode=xode, tode=tode, xr=xr, tr=tr)
    nothing
end


"Solve a given ODEProblem with Tsit5() algorithm"
function solve_ODE(prob::ODEProblem)
    return solve(prob, Tsit5())
end

"Solve a given SDEProblem with EM() algorithm and some dt"
function solve_SDE(sprob::SDEProblem; dt=1e-3) 
    return solve(sprob, EM(), dt=dt)
end

"Extract histograms on the distribution of the steady state"
function get_hist(sol::RODESolution; idx=1, tstart=2e2, bins=3e3:1e1:4.5e3)
    #~ Collect data
    # N = sum(sol.u[begin])
    x = reduce(hcat, sol.u)[idx,findfirst(sol.t .> tstart):end]
    fh = FHist.Hist1D(x; binedges=bins)
    return fh |> normalize
end

function store_hist(fh::FHist.Hist1D, fname::String)
    #~ Extract path from filename, and make sure it exists
    _path = rsplit(fname, "/", limit=2)
    mkpath(_path[begin])
    #~ Save 
    jldsave(fname; histogram = fh)
end

end # module SIRSReDynamics
#/ End module
