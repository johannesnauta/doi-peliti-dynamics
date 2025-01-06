#= Module for simulating dynamics of generalized (competitive) Lotka-Volterra model
 * includes functions that numerically integrate the reduced dynamics
 * if not specified, closed-form derivations of the reduced dynamics have been
   performed using Mathematica
=#
#/ Start module
module gLVReDynamics

#/ Packages
using Catalyst
using DifferentialEquations
using FHist
using JLD2
using Random
using Symbolics
using SparseArrays

using ModelingToolkit: t_nounits as t, D_nounits as D

include("../scripts/symbolics/LV.jl")
using .SymbolicLV

#################
### FUNCTIONS ###
"Solve generalized Lotka-Volterra for S species"
function solve_gLV(;
    S::Int = 24,
    μA::Float64 = 1.0,
    σA::Float64 = 0.1,
    r = ones(S),
    κ = 1e3 .* ones(S),
    c = 0.1,
    dt = 1e-2,
    tspan = (0.0,1e2),
    meanfield = false
    )
    sparserng = Random.Xoshiro(S*1994)
    spα = sprand(sparserng,S,S,c)
    rng = Random.Xoshiro(S*1234)
    #~ Rescale
    μA = μA / S
    σA = σA / sqrt(S)
    #~ sample interaction matrix (note: we do not care about the diagonal term)
    #~ note: in this simple case, we do not properly handle mutualistic interactions
    #        as such, just sample from a half-normal distribution
    α = abs.(μA .+ σA .* randn(rng,S,S)) .* spα

    rn = get_ReactionSystem(S)
    umap = [rn.x[i] => κ[i] ./ 10 for i in 1:S]
    pmap = [rn.r => r, rn.κ => κ, rn.α => α]
    sprob = SDEProblem(rn, umap, tspan, pmap)
    ssol = solve_SDE(sprob, dt=dt, seed=S*42)

    #/ Additionally solve mean-field eq. by solving the ODE
    if meanfield
        oprob = ODEProblem(rn, umap, tspan, pmap)
        osol = solve_ODE(oprob)
        return ssol, osol
    end
    return ssol
end

function get_ReactionSystem(S::Int)
    _t = default_t()
    @parameters r[1:S] κ[1:S]
    @parameters α[1:S,1:S]
    @species (x(_t))[1:S]
    rxns = []
    for i in 1:S
        push!(rxns, Reaction(r[i], [x[i]], [x[i]], [1], [2]))
        push!(rxns, Reaction(2*r[i]/κ[i], [x[i]], [x[i]], [2], [1]))
        for j in 1:S
            (j != i) && push!(rxns, Reaction(α[i,j]/κ[i], [x[i], x[j]], [x[j]], [1,1], [1]))
        end
    end
    @named rs = ReactionSystem(rxns, _t, [x...], [r..., κ..., α...])
    return complete(rs)
end

function solve_simple(S::Int; r=1e-1.*ones(S), κ=1e3.*ones(S), tspan=(0.0,5e2), noise=true)
    rn = get_ReactionSystem(S)
    umap = [rn.x[i] => κ[i] ./ 10 for i in 1:S]
    pmap = [rn.r => r, rn.κ => κ]
    if noise 
        sprob = SDEProblem(rn, umap, tspan, pmap)
        ssol = solve_SDE(sprob, dt=1e-2)
        return ssol
    end 
    prob = ODEProblem(rn, umap, tspan, pmap)
    sol = solve(prob, Tsit5(), save_everystep=false)
    return sol
end

########################
### REDUCED DYNAMICS ###
function solve_rgLV(
    drift, noise;
    S::Int = 24,
    focal::Int = 1,
    μA::Float64 = 1.0,
    σA::Float64 = 0.1,
    rv = ones(S),
    κv = 1e3 .* ones(S),
    tspan=(0.0,5e4)
    )
    sprob = get_rSDEProblem(
        drift, noise;
        S=S, focal=focal, μA=μA, σA=σA, rv=rv, κv=κv, tspan=tspan
    )
    ssol = solve_SDE(sprob, dt=1e-2, seed=S*54321)
    return ssol
end

function get_rSDEProblem(
    drift, noise
    ;
    S::Int = 24,
    focal::Int = 1,
    μA::Float64 = 1.0,
    σA::Float64 = 0.1,
    rv = ones(S),
    κv = 1e3 .* ones(S),
    tspan=(0.0,1e3)
    )
    rng = Random.Xoshiro(S*1234)
    #~ Rescale
    μA = μA / S
    σA = σA / sqrt(S)
    #~ sample interaction matrix (note: we do not care about the diagonal term)
    #~ note: in this simple case, we do not properly handle mutualistic interactions
    #        as such, just sample from a half-normal distribution
    αv = abs.(μA .+ σA .* randn(rng,S,S))
    #/ Define reduced SDEProblem
    @parameters r[1:S] κ[1:S] α[1:S,1:S]
    @variables η(t)
    deq = [D(η) ~ drift]
    seq = [sqrt(noise)]
    # @info "eqs" deq seq
    umap = [η => κv[focal]/10.0]
    pmap = [r => rv, κ => κv, α => αv]
    
    #/ Define ODESystem and SDESystem
    @named odesys = ODESystem(deq, t)
    @named sdesys = SDESystem(odesys, seq)
    sprob = SDEProblem(complete(sdesys), umap, tspan, pmap)
    # sprob = ODEProblem(complete(odesys), umap, tspan, pmap)
    return sprob
end

### Some simple functions for testing
"Define super simple single-species logistic growth"
function get_simpleReactionSystem()
    rn = @reaction_network begin
        r, X --> 2X
        2*r/κ, 2X --> X
    end
    return rn
end

"Solve super simple model"
function solve_simple(; r=1e-1, κ=1e3, tspan=(0.0,5e2), noise=true)
    rn = get_simpleReactionSystem()
    umap = [:X => 1e2]
    pmap = [:r => r, :κ => κ]
    if noise 
        sprob = SDEProblem(rn, umap, tspan, pmap)
        ssol = solve_SDE(sprob, dt=1e-2)
        return ssol 
    end 
    prob = ODEProblem(rn, umap, tspan, pmap)
    sol = solve(prob, Tsit5())
    return sol
end

########################
### HELPER FUNCTIONS ###
"Solve trajectories for the full system (stochastic & mean-field) and reduced system"
function solve_trajectories(; S=24, tspan=(0.0, 1e2), dt=1e-2, ddir="../data/trajectories/")
    ssol, osol = solve_gLV(; S=S, tspan=tspan, dt=dt, meanfield=true)
    drift, noise = SymbolicLV.get_driftnoise(S, focal=1)
    rsol = solve_rgLV(drift, noise; S=S, tspan=tspan)
    #/ Extract trajectories for i=1
    xsde = reduce(hcat, ssol.u)[1,:]
    tsde = ssol.t
    xode = reduce(hcat, osol.u)[1,:]
    tode = osol.t
    xr   = reduce(hcat, rsol.u)[1,:]
    tr   = rsol.t
    #~ Save
    fname = ddir*"gLVtrajectories.jld2"
    jldsave(fname; xsde=xsde, tsde=tsde, xode=xode, tode=tode, xr=xr, tr=tr)
    nothing
end

"Solve a given ODEProblem with Tsit5() algorithm"
function solve_ODE(prob::ODEProblem)
    return solve(prob, Tsit5())
end
    
"Solve a given SDEProblem with EM() algorithm and some dt"
function solve_SDE(sprob::SDEProblem; dt=1e-2, seed=42) 
    return solve(sprob, EM(), dt=dt, seed=seed)
end

"Extract histograms on the distribution of the steady state"
function get_hist(sol::RODESolution; idx=1, tstart=1e4, bins=7e2:5e0:1.1e3)
    #~ Collect data
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

end # module gLVReDynamics
#/ End module
