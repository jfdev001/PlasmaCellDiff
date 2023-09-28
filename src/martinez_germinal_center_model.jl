# Kinetic model of the gene regulatory module controlling the germinal 
# center (gc) exit pathway by using ODEs. 
# Implementation based on 
# Martínez MR, et. al. Quantitative modeling of the terminal differentiation of 
# B cells and mechanisms of lymphomagenesis. Proc Natl Acad Sci U S A. 2012 Feb 
# 14; 109(7):2672-7. (aka Martinez2012)

using DynamicalSystems: SVector, SMatrix, CoupledODEs
using UnPack

"""
    germinal_center_exit_pathway_rule(u, params::GerminalCenterODEParams, t)

Return rule for gene regulatory module controlling germinal center exit pathway
with coupled BCR and CD40 regulatory signals.

# Arguments 
- `u`: State vector of dynamical system. The elements of this vector represent
    protein levels of `p` (BLIMP1), `b` (BCL6), and `r` (IRF4).
- `params::GerminalCenterODEParams`: Typed mutable struct whose parametrized
    symbols are in `[:gaussian, :constant, :reciprocal]` and which determine
    the BCR/CD40 regulatory signaling mechanism.
- `t`: Current time for numerical integration. No need to pass an argument here
    since numerical integration is handled with builtin solvers.

# References
[1] : Equations S1 - S3 from Martinez2012
"""
function germinal_center_exit_pathway_rule(
    u, params::GerminalCenterODEParams, t)
    # parameters 
    @unpack μp, μb, μr = params
    @unpack σp, σb, σr = params
    @unpack kp, kb, kr = params
    @unpack λp, λb, λr = params

    # transcription factor state variables 
    p, b, r = u

    # compute scaled dissociation constants and protein levels
    kp_scaled = dissociation_scaler(kp, p)
    kb_scaled = dissociation_scaler(kb, b)
    kr_scaled = dissociation_scaler(kr, r)
    r_scaled = transcription_factor_scaler(kr, r)

    # regulatory signals
    bcr = BCR(u, params, t)
    cd40 = CD40(u, params, t)

    # system describing evolution of transcription factors in germinal center
    pdot = μp + σp*kb_scaled + σp*r_scaled - λp*p
    bdot = μb + σb*kp_scaled*kb_scaled*kr_scaled - (λb + bcr)*b
    rdot = μr + σr*r_scaled + cd40 - λr*r

    return SVector(pdot, bdot, rdot)
end

"""
    germinal_center_exit_pathway(u0, params::GerminalCenterODEParams) 

Return CoupledODEs model for gene regulatory module controlling germinal center 
exit pathway with coupled BCR and CD40 regulatory signals. 
"""
function germinal_center_exit_pathway(u0, params::GerminalCenterODEParams) 
    return CoupledODEs(germinal_center_exit_pathway_rule, u0, params)
end


"""
    bcr_subnetwork_rule(u, params::GerminalCenterODEParams, t)

BCR signaling rule decoupled from CD40 signaling and with negligible IRF4 level.

# References
[1] : Equations S6 and S7 from Martinez2012
"""
function bcr_subnetwork_rule(u, params::GerminalCenterODEParams, t)
    # parameters 
    @unpack μp, μb = params
    @unpack σp, σb = params
    @unpack kp, kb = params
    @unpack λp, λb = params
    @unpack bcr0 = params

    # BLIMP1 and BCL6 protein levels
    p, b = u

    # scaled dissociation constants
    kb_scaled = dissociation_scaler(kb, b)
    kp_scaled = dissociation_scaler(kp, p)

    # regulatory signal
    bcr = BCR(u, params, t)
    
    # BCR subnetwork ODEs
    pdot = μp + σp*kb_scaled - λp*p
    bdot = μb + σb*kp_scaled*kb_scaled - (λb + bcr)*b

    return SVector(pdot, bdot)
end 

function bcr_subnetwork(u0, params::GerminalCenterODEParams)
    return CoupledODEs(bcr_subnetwork_rule, u0, params)
end 
