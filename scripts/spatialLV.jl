#= Module for simulating (reduced) dynamics of spatial Lotka-Volterra model 
=#
#/ Start module
module sLVReDynamics

#/ Packages
using Catalyst
using DifferentialEquations
using Random
using Symbolics
using JLD2

# using ModelingToolkit: t_nounits as t, D_nounits as D

#################
### FUNCTIONS ###
function get_odeproblem(;
    N = 1.0,
    d = 0.2,
    c = 1.0,
    k = 0.8,
    ε = 0.3,
    u0 = N .* [(N - d)*0.25, (N-d)*0.25, (N-d)*0.5],
    tspan = (0.0, .5e2)
)
    rn = get_reactionsystem()    
    umap = [rn.A => u0[1], rn.B => u0[2], rn.E => u0[3]]
    pmap = [rn.N => N, rn.c => c, rn.k => k, rn.ε => ε]
    prob = ODEProblem(rn, umap, tspan, pmap)
    return prob
end

function get_sdeproblem(;
    N = 10000.0,
    d = 0.2,
    c = 1.0,
    k = 0.8,
    ε = 0.3,
    u0 = N.*[(1-d-0.05),0.05,0.05],
    tspan = (0.0, 1e2)
)
    rn = get_reactionsystem()
    umap = [rn.A => u0[1], rn.B => u0[2], rn.E => u0[3]]
    pmap = [rn.N => N, rn.c => c, rn.k => k, rn.ε => ε]
    sprob = SDEProblem(rn, umap, tspan, pmap)
    return sprob
end

function get_reactionsystem()
    _t = default_t()
    @parameters c k ε N
    @species A(_t) B(_t) E(_t)
    rxns = [
        Reaction(c/N, [A,E], [A], [1,1], [2]),
        Reaction(k/N, [A,B], [B], [1,1], [2]),
        Reaction(ε, [A], [E], [1], [1]),
        Reaction(ε, [B], [E], [1], [1])
    ]
    @named rs = ReactionSystem(rxns, _t, [A,B,E], [c,k,ε,N])
    return complete(rs)
end

function solve_odeproblem(prob::ODEProblem)
    return solve(prob, Tsit5())
end

function solve_sdeproblem(prob::SDEProblem; dt=1e-2, seed=1234)
    return solve(prob, EM(), dt=dt, seed=seed)
end

########################
### HELPER FUNCTIONS ###
"""Simple function to run and save dynamics for the three regions of interest
(1) Both species persist
(2) Only the prey species persists
(3) Both species go extinct
"""
function run_dynamics(; Dv=[0.1, 0.5, 0.9], ddir="../data/trajectories/")
    for i in eachindex(Dv)
        prob = get_odeproblem(; d=Dv[i])
        sol = solve_odeproblem(prob)
        #/ Extract states (only for A and B) and time
        x = reduce(hcat, sol.u)[1:2,:]
        t = sol.t
        #/ Save
        fname = ddir * "sLVtrajectory_D$(Dv[i]).jld2"
        jldsave(fname; x = x, t = t)
    end
    nothing
end



end # module sLVReDynamics
#/ End module
