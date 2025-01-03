#= Module for deriving Hamilton equations for the LV model =#
#/ Start module
module SymbolicLV

#/ Packages
using Groebner
using Symbolics
using SparseArrays
using Random 

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function get_driftnoise(S::Int; focal::Int=1, c=min(1.0,0.1))
    rng = Random.Xoshiro(S*1994)
    basesparr = sprand(rng,S,S,c)
    spα = similar(basesparr, Num)
    #~ Define variables 
    @parameters r[1:S] κ[1:S]
    @parameters α[1:S,1:S]
    @variables x[1:S] q[1:S]
    rows, cols, _ = findnz(basesparr)
    for (i,j) in zip(rows,cols)
        spα[i,j] = α[i,j]
    end

    #~ Construct Hamiltonian
    H = sum(r[i]*x[i]*(x[i]-1)*q[i] for i in 1:S) +
        sum(r[i]*x[i]*(1-x[i])*q[i]^2/κ[i] for i in 1:S) +
        sum(spα[i,j]/κ[i]*x[j]*(1-x[i])*q[i]*q[j] for i in 1:S, j in 1:S if j!=i)

    dn = Array{Equation}(undef, S)
    for j in 1:S
        dx = Differential(x[j])
        _dn = substitute(expand_derivatives(dx(H)), Dict([x[j] => 1 for j in 1:S]))
        #~ note: we can divide by q[j] by implicitly assuming that q[j]≠0
        #~       this makes the results only hold in the feasible regime
        #~       otherwise the solution cannot be solved easily, as it'll be non-linear
        dn[j] = simplify(_dn / q[j]) ~ 0
    end
    eqs = dn[filter(j->j!=focal, 1:S)]
    vars = [q[j] for j in 1:S if j!=focal]
    qsol = Symbolics.symbolic_linear_solve(eqs, vars)

    xdict = Dict(x[j] => 1 for j in 1:S if j!=focal)
    qdict = Dict(vars[i] => qsol[i] for i in eachindex(vars))
    Hred = substitute(H, merge(xdict, qdict))

    #~ Cole-Hopf variables
    @variables η(t) θ
    Hch = substitute(Hred, Dict(x[1] => exp(θ), q[1] => η*exp(-θ)))
    dθ1 = Differential(θ)^1
    dθ2 = Differential(θ)^2
    drift = substitute(expand_derivatives(dθ1(Hch)), Dict(θ => 0))
    noise = substitute(expand_derivatives(dθ2(Hch)), Dict(θ => 0))

    deq = [D(η) ~ drift]
    neq = [sqrt(noise)]
    # @named odesys = ODESystem(deq, t)
    
    return drift, noise
end


end # module SymbolicLV
#/ End module
