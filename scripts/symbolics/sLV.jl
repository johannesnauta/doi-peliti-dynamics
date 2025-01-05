#= Module for deriving Hamilton equations for the spatial LV model =#
#/ Start module
module SymbolicSpatialLV

#/ Packages
using Groebner
using Symbolics
using SparseArrays
using Random 

# using ModelingToolkit
# using ModelingToolkit: t_nounits as t, D_nounits as D

function get_reduced_Hamiltonian(; focal::Int=1)
    #~ Define variables 
    @parameters c k ε d N
    @variables x[1:3] q[1:3]

    #~ Construct Hamiltonian
    H = c / N * x[1] * (x[1] - x[3]) * q[1] * q[3] + ε * (x[3] - x[1]) * q[1] +
        k / N * x[2] * (x[2] - x[1]) * q[1] * q[2] + ε * (x[3] - x[2]) * q[2]

    dn = Array{Equation}(undef, length(x) + 1)
    for i in eachindex(x)
        dx = Differential(x[i])
        deriv = expand_derivatives(dx(H))
        subdict = Dict([x[i] => 1 for i in 1:length(x)])
        _dn = simplify(substitute(deriv, subdict))
        dn[i] = _dn ~ 0
    end
    dn[end] = N - d - q[1] - q[2] - q[3] ~ 0
    @info "mean-field eqs" dn
        
    eqs = dn[filter(j->j!=focal, 1:length(dn))]
    vars = [q[j] for j in 1:length(x) if j!=focal]
    # ~ There are two solution, one where qB==0, and one where it's not
    eqs = vcat(dn[3], dn[4])
    vars = [q[2], q[3]]
    
    qsol = simplify.(Symbolics.symbolic_linear_solve(eqs, vars))
    @info "solution" qsol

    xdict = Dict(x[j] => 1 for j in 1:length(x) if j!=focal)
    qdict = Dict(vars[i] => qsol[i] for i in eachindex(vars))
    Hred = substitute(H, merge(xdict, qdict))

    #~ Cole-Hopf variables
    @variables η θ
    Hch = substitute(Hred, Dict(x[1] => exp(θ), q[1] => η*exp(-θ)))    
    return Hch
end

# function get_field(H)    
#     @variables η θ
#     ∂η = Differential(η)
#     ∂θ = Differential(θ)
#     dη = simplify(expand_derivatives(∂θ(H)))
#     dθ = simplify(expand_derivatives(-∂η(H)))
#     return dη, dθ
# end

# function get_driftnoise(H)
#     @variables θ
#     dθ1 = Differential(θ)^1
#     dθ2 = Differential(θ)^2
#     drift = substitute(expand_derivatives(dθ1(H)), Dict(θ => 0))
#     noise = substitute(expand_derivatives(dθ2(H)), Dict(θ => 0))
#     return drift, noise
# end

function get_field()
    @variables η θ
    @variables ε c k D
    Hstar = exp(-θ) * (exp(θ) - 1) * η * (ε - c * ((1 - D)*exp(θ) - η))
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = simplify(expand_derivatives(∂θ(Hstar)))
    dθ = simplify(expand_derivatives(-∂η(Hstar)))
    return dη, dθ
end



end # module SymbolicLV
#/ End module
