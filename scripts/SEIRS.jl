#= Module for simulating the SEIRS model with two E states =#
#/ Start module
module SEIRS

#/ Packages
using Catalyst
using DifferentialEquations
using FHist

#################
### FUNCTIONS ###
function define_JumpProblem(
    ;
    Nval::Int = 10_000,    
    U0::Array = ones(Int, 5),
    βval::Float64 = 1.1094736842,
    μval::Float64 = 0.2,
    vval::Float64 = 0.4,
    wval::Float64 = 0.1,
    pval::Float64 = 0.9,
    γval::Float64 = 0.3,
    tspan = (0.0, 10_000.),
    kwargs...
    )
	  @variables t
    @species S(t) E1(t) E2(t) I(t) R(t)
    @parameters β γ μ p v w N

    #/ Specify reactions
    rxns = [
        Reaction(μ, [S], nothing),
        Reaction(μ, [E1], nothing),
        Reaction(μ, [E2], nothing),
        Reaction(μ, [I], nothing),
        Reaction(μ, [R], nothing),
        Reaction(μ*N, nothing, [S]),
        Reaction(p*β/N, [S, I], [E1, I], [1,1], [1,1]),
        Reaction((1-p)*β/N, [S, I], [E2, I], [1,1], [1,1]),
        Reaction(v, [E1], [I]),
        Reaction(w, [E2], [I]),
        Reaction(γ, [I], [R])
    ]
    #/ Create reaction system
    @named rs = ReactionSystem(rxns, t)

    #/ Specify parameters
    params = [N => Nval, β => βval, μ => μval, v => vval, w => wval, p => pval, γ => γval]
    u0 = [S => U0[1], E1 => U0[2], E2 => U0[3], I => U0[4], R => U0[5]]
    
    #/ Define DiscreteProblem
    dprob = DiscreteProblem(rs, u0, tspan, params)
    jprob = JumpProblem(rs, dprob, Direct(), save_positions=(false,false))
    return jprob
end

"Compute the histogram of the steady state abundances of the susceptible population"
function get_Shist(
    ;
    ntrajectories = 50,
    tspan = (0.0, 5_000.0),
    kwargs...
    )
	  #/ Get EnsembleProblem
    jprob = define_JumpProblem(; tspan=tspan, kwargs...)
    output_func(sol, i) = (sol[1,:], false)
    eprob = EnsembleProblem(jprob, output_func=output_func)
    #/ Define times at which to save state
    #~ Note: can also do that within the output func, but do it like this for now
    tstart = tspan[end] / 4 * 3
    tpoints = collect(LinRange(tstart, tspan[end], 100))
    #/ Solve
    esol = solve(
        eprob, SSAStepper(), EnsembleThreads(); trajectories=ntrajectories, saveat=tpoints
    )
    #/ Collect trajectories and generate histogram
    _esols = [sol[begin+1:end] for sol in esol]
    states = reduce(vcat, _esols)
    fh = FHist.Hist1D(states, 6000:25:8000)
    return fh |> FHist.normalize
end

function define_ODEProblem(
    ;
    N::Int = 10_000,    
    Q0::Array = ones(Int, 5),
    βval::Float64 = 0.789,
    μval::Float64 = 0.2,
    vval::Float64 = 0.4,
    wval::Float64 = 0.1,
    pval::Float64 = 0.9,
    γval::Float64 = 0.3,
    tspan = (0.0, 10_000.)
    )
    @variables t
    @variables (q(t))[1:5]
    @parameters β γ μ p v w
    D = Differential(t)
    
    # Notes: v = ν1, w = ν2
	  eqs = [
        D(q[1]) ~ -(1-p)*q[1]*q[4]*β/N - p*q[1]*q[4]*β/N + N*μ - q[1]*μ,
        D(q[2]) ~ p*q[1]*q[4]*β/N - q[2]*μ - q[2]*v,
        D(q[3]) ~ (1-p)*q[1]*q[4]*β/N - q[3]*μ - q[3]*w,
        D(q[4]) ~ -q[4]*γ - q[4]*μ + q[2]*v + q[3]*w,
        D(q[5]) ~ q[4]*γ - q[5]*μ
    ]
    @named sys = ODESystem(eqs, t)

    pmap = [β => βval, μ => μval, v => vval, w => wval, p => pval, γ => γval]
    qmap = [q[i] => Q0[i] for i in 1:5]

    prob = ODEProblem(sys, qmap, tspan, pmap)
    return prob
end

end # module SEIRS
#/ End module
