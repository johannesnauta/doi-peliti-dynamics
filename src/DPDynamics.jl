#/ Start module
module DPDynamics

#/ Export

#/ Macro to make the LSP (Eglot) aware of this module within some scripts
#  that would otherwise not be properly linted.
#~ Other files wherein the LSP should recognize the HostPathogen module can
#~ be included here.
#!Note: this may safely be ignored when not using Emacs+Eglot
macro ignore(args...) end
#/ 'include' files by ignoring
@ignore include("../scripts/reduced-dynamics.jl")
@ignore include("../plot/plot.jl")

end # module DPDynamics
#/ End module
