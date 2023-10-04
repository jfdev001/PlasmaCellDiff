# `ChaosTools.fixedpoints` uses `t=0.0` during root finding. I defined
# functions here that override the existing root finding helper functions
# so that a non-zero timestep invokes gaussian signaling of BCR/CD40
using PlasmaCellDiff
import ChaosTools: to_root_f, to_root_J

function to_root_f(
    ds::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing)
    return u -> dynamic_rule(ds)(u, p, p.regulation_timestep_t)
end 

function to_root_J(
    Jf, ::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing) 
    return u -> Jf(u, p, p.regulation_timestep_t)
end
