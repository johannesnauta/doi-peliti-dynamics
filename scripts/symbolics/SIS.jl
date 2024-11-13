#= Module for deriving Hamilton equations for the SIS model =#
#/ Start module
module SymbolicSIS

#/ Packages
using Symbolics

#################
### FUNCTIONS ###
function getfield(; modified=false)
    (modified) && (return getmodifiedfield())
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

function getmodifiedfield()
    #/ Define Hamiltonian
    @variables η θ
    @variables N γ ρ α
    H = (1 - exp(θ)) * (N - η) * (α*η*exp(-θ)*(N-η)/N^2 + γ*η*exp(-θ)/N - ρ)
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = simplify(expand_derivatives(∂θ(H)))
    dθ = simplify(expand_derivatives(-∂η(H)))
    return (dη, dθ)
end

end # module SymbolicSIS
#/ End module
