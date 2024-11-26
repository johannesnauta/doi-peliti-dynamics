#= Module for functions that solve and/or integrate Levins model =#
#/ Start module
module Levins

#/ Packages
using DifferentialEquations, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

#################
### FUNCTIONS ###
function generate_fixed_ODESystem()
    @parameters c[1:2] e[1:2] d
    @variables (n(t))[1:3]
    eqns = [
        D(n[1]) ~ c[1]*n[1]*(1 - d - n[1] - n[2]) - e[1]*n[1],
        D(n[2]) ~ c[2]*n[1]*n[2] - e[2]*n[2],
        D(n[3]) ~ -c[1]*n[1]*(1-d-n[1]-n[2]) - c[2]*n[1]*n[2] + e[1]*n[1] + e[2]*n[2]
    ]
    @named sys = ODESystem(eqns, t)
    return complete(sys)
end

function generate_ODESystem(; S::Int = 2)
    @parameters c[1:S] e[1:S] d
    @variables (n(t))[1:(S+1)]
    eqns = [D(n[i]) ~ c[i]*n[i]*(n[end] - d) - e[i]*n[i] for i in 1:S]
    eqn  = [D(n[end]) ~ -sum(eq.rhs for eq in eqns)]
    eqns = [eqns..., eqn...]
    @named sys = ODESystem(eqns, t)
    return complete(sys)
end

function generate_ODEProblem(sys;
    c=[1.0,0.5],
    e=[0.2,0.2],
    d=0.0,
    u=[0.99, 0.01, 0.0],
    tspan=(0.0, 1e3)
    )
    @info "npred" 1/c[1] * (1 - d - e[1] - e[2]/c[2]) 1-d-e[1]/c[1] 1-e[1]/c[1]
    pmap = [sys.c => c, sys.e => e, sys.d => d]
    xmap = [sys.n => u]
    prob = ODEProblem(sys, xmap, tspan, pmap)
    return prob
end

end # module Levins
#/ End module
