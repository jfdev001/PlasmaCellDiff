# Functions for computing nullclines of equations defined in Martinez2012
# to model germinal cell center regulation

using PlasmaCellDiff: GerminalCenterODEParams, BCR, CD40
using PlasmaCellDiff: dissociation_scaler, transcription_factor_scaler
using UnPack

function blimp1_steady_state(u, params::GerminalCenterODEParams)     
    # parameters 
    @unpack μp, μb, μr = params
    @unpack σp, σb, σr = params
    @unpack kp, kb, kr = params
    @unpack λp, λb, λr = params

    # transcription factor state variables 
    _, b, r = u

    # pdot = 0 = solve[pdot, p]
    kb_scaled = dissociation_scaler(kb, b) 
    r_scaled = transcription_factor_scaler(kr, r) 
    p_steady_state = (μp + σp*kb_scaled + σp*r_scaled)/(λp)
    
    return p_steady_state 
end 

function bcl6_steady_state(u, params::GerminalCenterODEParams) 
    # parameters 
    @unpack μp, μb, μr = params
    @unpack σp, σb, σr = params
    @unpack kp, kb, kr = params
    @unpack λp, λb, λr = params

    # transcription factor state variables 
    p, b, r = u

end 

function irf4_steady_state(u, params::GerminalCenterODEParams)
    # parameters 
    @unpack μp, μb, μr = params
    @unpack σp, σb, σr = params
    @unpack kp, kb, kr = params
    @unpack λp, λb, λr = params

    # transcription factor state variables 
    p, b, r = u
end 
