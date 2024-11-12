#= Module for deriving Hamilton equations for the SIS model =#
#/ Start module
module SymbolicSIS

#/ Packages
using Symbolics

#################
### FUNCTIONS ###
function getfield()
    #/ Define Hamiltonian
    @variables η θ
    @variables N γ ρ
    H = (1 - exp(θ)) * (N - η) * ((γ*η/N)*exp(-θ) - ρ)
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = simplify(expand_derivatives(∂θ(H)))
    dθ = simplify(expand_derivatives(-∂η(H)))
    return dη, dθ
end

end # module SymbolicSIS
#/ End module
