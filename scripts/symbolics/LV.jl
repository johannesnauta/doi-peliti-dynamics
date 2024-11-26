#= Module for deriving Hamilton equations for the LV model =#
#/ Start module
module SymbolicLV

#/ Packages
using Symbolics 

function getfield()
    #/ Define Hamiltonian
    @variables η θ
    @variables N D c k ε
    # H = η * (1 - exp(θ)) * ( c * (N - D - η) - ε*exp(-θ))
    H = η * (1 - exp(θ)) * (ε*exp(-θ) - ε - k*(N-D-η-ε/c))
    # H = exp(-θ)*(1-exp(θ))*ε*η + (exp(θ) - 1)*ε*η + (exp(θ) - 1)*k*(N-D-ε/c-η)*η
    ∂η = Differential(η)
    ∂θ = Differential(θ)
    dη = simplify(expand_derivatives(∂θ(H)))
    dθ = simplify(expand_derivatives(-∂η(H)))
    return dη, dθ
end

end # module SymbolicLV
#/ End module
