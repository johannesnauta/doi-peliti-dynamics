#/ Start module
module Dynamics

#/ Packages
using Catalyst
using DifferentialEquations
using FHist
using Random

#################
### EPIDEMICS ###
function define_SIRSrs()
	  @variables t
    @species S(t) I(t) R(t)
    @parameters β γ α σ

    nreactions = 4
    rxs = Array{Reaction}(undef, nreactions)

    rxs[1] = Reaction(β/(S+I+R), [S, I], [I], [1, 1], [2])
    rxs[2] = Reaction(γ, [I], [R])
    rxs[3] = Reaction(α/(S+I+R)^2, [S, I], [I], [1, 2], [3])
    rxs[4] = Reaction(σ, [R], [S])
    @named rs = ReactionSystem(rxs)
    return rs
end

function get_SIRSparams(rs::ReactionSystem)
    #/ Get parameters
    @unpack β, γ, α, σ = rs
    params = (β => 4/5, α => 1/15, γ => 2/5, σ => 1)
    return params
end

function get_SIRSode(
    ;
    tspan = (0.0, 1000.0),
    X0 = [900., 100., 0.]
    )
    rs = define_SIRSrs()
    params = get_SIRSparams(rs)
    #/ Define ODE system
    odeprob = ODEProblem(rs, X0, tspan, params)
    return odeprob
end

function get_SIRSsde(
    ;
    tspan = (0.0, 100.0),
    X0 = [900., 100., 0.]
    )
    rs = define_SIRSrs()
    params = get_SIRSparams(rs)
    #/ Define ODE system
    sdeprob = SDEProblem(rs, X0, tspan, params)
    return sdeprob
end

function get_reduced_SIRSsde(; X0=[900., 100., 0.], tspan=(0.0, 1000.0))
    β = 4/5
    α = 1/15
    γ = 2/5
    σ = 1.0
    params = (β, α, γ, σ)
    n = sum(X0)
    
    function f(u, p, t)
        #~ Note: not generic as Mathematica...
        β, γ, α, σ = p
        return (
            -((5/7) * (-400 + (2 * u)/5)) +
                ((-1000 + u) * u) / 1750 -
                ((-1000 + u)^2 * u) / 29400000
        )        
    end
    function g(u, p, t)        
        #~ Note: not generic as Mathematica...
        return sqrt(
            -(5/7) * (-400 + (2 * u)/5) -
                ((-1000 + u) * u) / 1750 +
                ((-1000 + u)^2 * u) / 29400000
        )
    end
    
    red_prob = SDEProblem(f, g, X0[begin], tspan, params)
    return red_prob
end

### Get histograms
function SIRS_full_abundance_hist(focal::Int; tspan=(0.0, 10_000.0), dt=0.01)
    Random.seed!(42)
	  sdeprob = get_SIRSsde(tspan=tspan)
    sdesol = solve(sdeprob, saveat=dt)
    tstart = trunc(Int, 50 / dt)
    fsol = sdesol[focal,:][tstart:end]
    h = FHist.Hist1D(fsol, 355:10:645)
    return h |> FHist.normalize
end

function SIRS_reduced_abundance_hist(; tspan=(0.0, 10_000.0), dt=0.01)
	  Random.seed!(42)
    redprob = get_reduced_SIRSsde(tspan=tspan)
    sdesol = solve(redprob, saveat=dt)
    tstart = trunc(Int, 50 / dt)
    fsol = sdesol.u[tstart:end]
    h = FHist.Hist1D(fsol, 355:10:645)
    return h |> FHist.normalize
end

######################
### LOTKA-VOLTERRA ###
"Define Lotka-Volterra reaction system"
function define_LVrs(n::Int)
    @variables t
    @species (x(t))[1:n]
    @parameters r[1:n] λ[1:n] A[1:n,1:n]

    nreactions = 2*n + n^2
    rxs = Array{Reaction}(undef, nreactions)

    _idx = 1
    for i in 1:n
        rxs[_idx] = Reaction(r[i], [x[i]], [x[i]], [1], [2])
        _idx += 1
        rxs[_idx] = Reaction(λ[i], [], [x[i]])
        _idx += 1
        for j in 1:n
            if i == j
                rxs[_idx] = Reaction(2A[i,j], [x[i]], [x[j]], [2], [1])
            else
                rxs[_idx] = Reaction(A[i,j], [x[i],x[j]], [x[j]], [1,1], [1])
            end
            _idx += 1
        end
    end
    @named rs = ReactionSystem(rxs)
    return rs    
end

"Get specific set of parameters for Lotka-Volterra model"
function get_LVparams(rs::ReactionSystem)	  
    n = 8    
	  Aval = [
        1 1 1/2 0 3/10 0 0 0;
        1/15 11/10 0 0 1/17 0 0 0;
        0 0 8/10 0 0 1 0 0;
        0 0 0 1 0 0 0 1/9;
        0 0 0 0 11/10 0 0 0;
        1/10 0 0 0 0 1 0 0;
        0 0 0 0 0 0 9/10 0;
        5/30 0 0 0 0 0 2/10 1
    ]' ./ 1000
    rval = 1/2
    λval = 1/2

    #/ Get parameters
    @unpack r, λ, A = rs
    rparams = (r[i] => rval for i in 1:n)
    λparams = (λ[i] => λval for i in 1:n)
    Aparams = (A[i,j] => Aval[i,j] for i in 1:n, j in 1:n)
    params = (rparams..., λparams..., Aparams...)
    return params
end

### Get ODEProblems
function get_simple_reduced_ode(
    ;
    tspan = (0.0, 100.0),
    X0 = 1000 .* ones(8)
    )
	  f(u, p, t) = 1/10 + u/2 - u*u/100
    red_prob = ODEProblem(f, 100., tspan)
    return red_prob
end

### Get SDEProblems
function get_reduced_sde(
    ;
    tspan = (0.0, 50.0),
    X0 = 1000.
    )
    #/ Define reduced (deterministic) dynamics
    f(u,p,t) = 1/2 + u/2 -
        ((425 - 3*sqrt(70) + 5*sqrt(7315 - 102*sqrt(70))) * u) / 10800 -
        u^2 / 1000 -
        (u * (500 - u + sqrt(252200 - 1000 * u + u^2))) / 33000 -
        (u * (
            500 - 5/16 * (1000 - u + sqrt(1006400 - 2000 * u + u^2)) +
                sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000 * u + u^2)))^2)
        )) / 20000

    g(u,p,t) = sqrt(
        1/2 + u/2 +
        ((425 - 3*sqrt(70) + 5*sqrt(7315 - 102*sqrt(70))) * u) / 10800 +
        u^2 / 1000 -
        (u * (u + (1000*u - 2*u^2) / (2*sqrt(252200 - 1000*u + u^2)))) / 16500 +
        (u * (500 - u + sqrt(252200 - 1000*u + u^2))) / 33000 -
        (u * (
            -(5/16) * (u + (2000*u - 2*u^2) / (2*sqrt(1006400 - 2000*u + u^2))) +
                (5 * (u + (2000*u - 2*u^2) / (2*sqrt(1006400 - 2000*u + u^2))) *
                (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))) /
                (16 * sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))^2))
        )) / 10000 +
            (u * (
                500 - 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)) +
                    sqrt(2000 + (-500 + 5/16 * (1000 - u + sqrt(1006400 - 2000*u + u^2)))^2)
            )) / 20000
    )
    
    red_prob = SDEProblem(f, g, X0, tspan)
    return red_prob
end

############
# WRAPPERS #
"Get ODEProblem for Lotka-Volterra system with 8 species"
function get_LVode(
    ;
    tspan = (0.0, 100.0),
    X0 = 1000 .* ones(8)
    )
    rs = define_LVrs(length(X0))
    params = get_LVparams(rs)
    #/ Define ODE system
    odeprob = ODEProblem(rs, X0, tspan, params)
    return odeprob
end

"Get SDEProblem for Lotka-Volterra system with 8 species"
function get_LVsde(
    ;
    n::Int = 8,
    tspan = (0.0, 50.0),
    X0 = 1000 .* ones(n)
    )
    rs = define_LVrs(length(X0))
    params = get_LVparams(rs)
    #/ Define ODE system
    sdeprob = SDEProblem(rs, X0, tspan, params)
    return sdeprob
end

"Get histogram by solving the full SDE and binning states in the stationary regime"
function full_abundance_hist(focal::Int; tspan=(0.0, 1_000.0), dt=0.01)
    Random.seed!(42)
	  sdeprob = get_LVsde(tspan=tspan)
    sdesol = solve(sdeprob, saveat=dt)
    tstart = trunc(Int, 50 / dt)
    fsol = sdesol[focal,:][tstart:end]
    h = FHist.Hist1D(fsol, 300:5:500)
    return h |> FHist.normalize
end

"Get histogram by solving reduced SDE and collecting states in the stationary regime"
function reduced_abundance_hist(; tspan=(0.0, 10_000.0), dt=0.01)
	  Random.seed!(42)
    redprob = get_reduced_sde(tspan=tspan)
    sdesol = solve(redprob, saveat=dt)
    tstart = trunc(Int, 50 / dt)
    fsol = sdesol.u[tstart:end]
    h = FHist.Hist1D(fsol, 300:5:500)
    return h |> FHist.normalize
end

end # module redDynamics
#/ End module
