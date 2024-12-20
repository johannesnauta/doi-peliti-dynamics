#= Module for simulating (reduced) dynamics of spatial Lotka-Volterra model 
=#
#/ Start module
module sLVReDynamics

#/ Packages
using Catalyst
using DifferentialEquations
using Random
using Symbolics

using ModelingToolkit: t_nounits as t, D_nounits as D

#################
### FUNCTIONS ###
function get_sdeproblem(;
    N = 1.0,
    c = 1.0,
    k = 0.4,
    ε = 0.3,
    u0 = [0.2, 0.2, 0.6],
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
        Reaction(ε, [A], []),
        Reaction(ε, [B], [])
    ]
    @named rs = ReactionSystem(rxns, _t)
    return complete(rs)
end

function solve_sdeproblem(prob::SDEProblem; dt=1e-2, seed=1234)
    return solve(prob, EM(), dt=dt, seed=seed)
end



end # module sLVReDynamics
#/ End module
