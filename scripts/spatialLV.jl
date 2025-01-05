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
function get_odeproblem(;
    N = 1.0,
    d = 0.2,
    c = 1.0,
    k = 0.8,
    ε = 0.3,
    u0 = [N*(1-d-0.05),0.05,0.05],
    tspan = (0.0, 1e2)
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
    @info "init cond." u0
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
    @named rs = ReactionSystem(rxns, _t)
    return complete(rs)
end

function solve_odeproblem(prob::ODEProblem)
    return solve(prob, Tsit5())
end


function solve_sdeproblem(prob::SDEProblem; dt=1e-2, seed=1234)
    return solve(prob, EM(), dt=dt, seed=seed)
end



end # module sLVReDynamics
#/ End module
