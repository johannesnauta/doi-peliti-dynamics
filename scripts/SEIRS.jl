#= Module for simulating the SEIRS model with two E states =#
#/ Start module
module SEIRS

#/ Packages
using Catalyst
using DifferentialEquations

#################
### FUNCTIONS ###
function define_DiscreteProblem()
	  @variables t
    @species S(t) E1(t) E2(t) I(t) R(t)
    @parameters β γ μ p v w N

    #/ Birth
    rxns = [
        Reaction(μ*N, nothing, [S]),
        Reaction(μ, [S], nothing),
        Reaction(μ, [E1], nothing),
        Reaction(μ, [E2], nothing),
        Reaction(μ, [I], nothing),
        Reaction(μ, [R], nothing),
        Reaction(β/N, [S, I], [E1, I], [1,1], [1,1]),
        Reaction(β/N, [S, I], [E2, I], [1,1], [1,1]),
        Reaction(v, [E1], [I]),
        Reaction(w, [E2], [I]),
        Reaction(γ, [I], [R])
    ]
    @named rs = ReactionSystem(rxns, t)
    return rs
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

function define_JumpProblem(
    ;
    N::Int = 10_000,
    )
	  sys = define_system(N)
    pmap = [β => 0.1]
    return sys, pmap
end

end # module SEIRS
#/ End module
