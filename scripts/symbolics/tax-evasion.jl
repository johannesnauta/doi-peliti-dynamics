#= Module for deriving Hamilton equations for the tax-evasion model =#
#/ Start module
module SymbolicTE

#/ Packages
using Symbolics

#################
### FUNCTIONS ###
function getfield()
    #/ Define Hamiltonian
    @variables η θ
    @variables N λ α β δ ε 
    H = (η - 1) * (1 - exp(θ)) * (exp(θ)*(α*(β+δ*η) + β*ε*η) + η*(δ*ε*η - λ*α)) /
        ((α + β) + δ*η)
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = simplify(expand_derivatives(∂θ(H)))
    dθ = simplify(expand_derivatives(-∂η(H))) 
    return dη, dθ
end

end # module SymbolicTE
#/ End module
