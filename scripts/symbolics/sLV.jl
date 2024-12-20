#= Module for deriving Hamilton equations for the spatial LV model =#
#/ Start module
module SymbolicSpatialLV

#/ Packages
using Groebner
using Symbolics
using SparseArrays
using Random 

using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

function get_driftnoise(; focal::Int=1)
    #~ Define variables 
    @parameters c k ε D N
    @variables x[1:3] q[1:3]

    #~ Construct Hamiltonian
    H = c / N * x[1] * (x[1] - x[3]) * q[1] * q[3] + ε * (x[3] - x[1]) * q[1] +
        k / N * x[2] * (x[2] - x[1]) * q[1] * q[2] + ε * (x[3] - x[2]) * q[2]

    dn = Array{Equation}(undef, length(x) + 1)
    for i in eachindex(x)
        dx = Differential(x[i])
        deriv = expand_derivatives(dx(H))
        subdict = Dict([x[i] => 1 for i in 1:length(x) if i != focal])
        _dn = simplify(substitute(deriv, subdict))
        dn[i] = _dn ~ 0
    end
    dn[end] = N - D - q[1] - q[2] - q[3] ~ 0
    @info "mean-field eqs" dn
        
    # eqs = dn[filter(j->j!=focal, 1:length(dn))]
    # vars = [q[j] for j in 1:length(x) if j!=focal]
    eqs = vcat(dn[3], dn[4])
    vars = [q[2], q[3]]
    
    qsol = simplify.(Symbolics.symbolic_linear_solve(eqs, vars))
    @info "solution" qsol

    xdict = Dict(x[j] => 1 for j in 1:length(x) if j!=focal)
    qdict = Dict(vars[i] => qsol[i] for i in eachindex(vars))
    Hred = substitute(H, merge(xdict, qdict))
    return Hred

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
