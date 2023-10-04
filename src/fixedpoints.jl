# `ChaosTools.fixedpoints` uses `t=0.0` during root finding. I defined
# functions here that override the existing root finding helper functions
# so that a non-zero timestep invokes gaussian signaling of BCR/CD40
using PlasmaCellDiff
import ChaosTools

function to_root_f(
    ds::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing)
    return u -> dynamic_rule(ds)(u, p, p.bcr0_t)
end 

function to_root_J(
    Jf, ::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing) 
    return u -> Jf(u, p, p.cd0_t)
end 
