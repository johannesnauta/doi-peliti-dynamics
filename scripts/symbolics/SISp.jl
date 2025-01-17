#= Module for deriving Hamilton equations for the SIS+ model, which is a modified SIS model
   with an additional reaction S+2I→3I with rate α
=#
#/ Start module
module SymbolicSIS

#/ Packages
using Symbolics

#################
### FUNCTIONS ###
function getfield()
    #/ Define Hamiltonian
    @variables η θ
    @variables N γ ρ α
    H = (1 - exp(θ)) * (N - η) * (η * exp(-θ) * (N*γ + α*(N-η)) / N^2 - ρ)
        # (α*η*exp(-θ)*(N-η)/N^2 + γ*η*exp(-θ)/N - ρ)
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = expand_derivatives(∂θ(H))
    dθ = expand_derivatives(-∂η(H))
    return (dη, dθ)
end

end # module SymbolicSIS
#/ End module
