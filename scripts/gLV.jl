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

#################
### FUNCTIONS ###
function solve_gLV(;
    S::Int = 32,
    μA::Float64 = 2.0,
    σA::Float64 = 0.1,
    r = ones(S),
    κ = 1e3 .* ones(S),
    tspan=(0.0,5e2)
    )
    rng = Random.Xoshiro(S*1234)
    #~ Rescale
    μA = μA / S
    σA = σA / sqrt(S)
    #~ sample interaction matrix (note: we do not care about the diagonal term)
    #~ note: in this simple case, we do not properly handle mutualistic interactions
    #        as such, just sample from a half-normal distribution
    α = abs.(μA .+ σA .* randn(rng,S,S))

    rn = get_ReactionSystem(S)
    umap = [rn.x[i] => κ[i] ./ 10 for i in 1:S]
    pmap = [rn.r => r, rn.κ => κ, rn.α => α]
    sprob = SDEProblem(rn, umap, tspan, pmap)
    ssol = solve_SDE(sprob, dt=1e-2, seed=S*2345)
    return ssol
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

function get_ReactionSystem(S::Int)
    t = default_t()
    @parameters r[1:S] κ[1:S]
    @parameters α[1:S,1:S]
    @species (x(t))[1:S]
    rxns = []
    for i in 1:S
        push!(rxns, Reaction(r[i], [x[i]], [x[i]], [1], [2]))
        push!(rxns, Reaction(2*r[i]/κ[i], [x[i]], [x[i]], [2], [1]))
        for j in 1:S
            (j != i) && push!(rxns, Reaction(α[i,j]/κ[i], [x[i], x[j]], [x[j]], [1,1], [1]))
        end
    end
    @named rs = ReactionSystem(rxns, t)
    return complete(rs)
end


function get_simpleReactionSystem()
    rn = @reaction_network begin
        r, X --> 2X
        2*r/κ, 2X --> X
    end
    return rn
end

"Solve a given SDEProblem with EM() algorithm and some dt"
function solve_SDE(sprob::SDEProblem; dt=1e-3, seed=42) 
    return solve(sprob, EM(), dt=dt, saveat=10.0, seed=seed)
end

end # module gLVReDynamics
#/ End module
