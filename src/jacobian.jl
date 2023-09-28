using UnPack
include("common.jl")

"""
    germinal_center_exit_pathway_jacobian(
        u, params::GerminalCenterODEParams{:reciprocal, :reciprocal}, t)

Return Jacobian of `germinal_center_exit_pathway_rule` w/ reciprocal BCR/CD40.

See `notebooks/germinal_cell_jacobian.pdf` for the source that generated the
symbolic Jacobian.

TODO: Check if both are requried to be reciprocal for this to be true.
"""
function germinal_center_exit_pathway_jacobian(
    u, params::GerminalCenterODEParams{:reciprocal, :reciprocal}, t)
    @unpack μp, μb, μr = params 
    @unpack σp, σb, σr = params
    @unpack kp, kb, kr = params 
    @unpack λp, λb, λr = params 
    @unpack bcr0, cd0, C0 = params

    p, b, r = u

    J_11 = -λp
    J_12 = -(2*b*kb^2*σp)/(b^2 + kb^2)^2 
    J_13 = -(2*r^3*σp)/(r^2 + kr^2)^2 + (2*r*σp)/(r^2 + kr^2) 

    J_21 = -(2*p*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)^2*(r^2 + kr^2))
    J_22 = (2*b^2*bcr0*kb^2)/(b^2+kb^2)^2 - (bcr0*kb^2)/(b^2 + kb^2) - 
        λb - (2*b*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)^2*(p^2 + kp^2)*(r^2 + kr^2))
    J_23 = -(2*r*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)*(r^2 + kr^2)^2)
 
    J_31 = 0
    J_32 = -(2*b*cd0*kb^2)/(b^2 + kb^2)^2
    J_33 = -λr - (2*r^3*σr)/(r^2 + kr^2)^2 + (2*r*σr)/(r^2 + kr^2)

    J = [J_11 J_12 J_13
         J_21 J_22 J_23
         J_31 J_32 J_33]

    return SMatrix{3, 3}(J)
end 

"""
    germinal_center_exit_pathway_jacobian(
        u, params::GerminalCenterODEParams{:gaussian, :gaussian}, t) 

Return Jacobian of `germinal_center_exit_pathway_rule` w/ gaussian BCR/CD40.

See `notebooks/germinal_cell_gaussian_jacobian.pdf` for the source that 
generated the symbolic Jacobian.

TODO: check if both are required to be gaussian for this Jacobian to be correct.
"""
function germinal_center_exit_pathway_jacobian(
    u, params::GerminalCenterODEParams{:gaussian, :gaussian}, t) 
    # parameters 
    @unpack μp, μb, μr = params 
    @unpack σp, σb, σr = params 
    @unpack kp, kb, kr = params 
    @unpack λp, λb, λr = params 
    @unpack bcr0, cd0, C0 = params
   
    # gaussian regulation parameters 
    @unpack bcr_max_signal = p
    @unpack bcr_max_signal_centered_on_timestep = p
    @unpack bcr_max_signal_timestep_std = p

    @unpack cd40_max_signal = p
    @unpack cd40_max_signal_centered_on_timestep = p
    @unpack cd40_max_signal_timestep_std = p

    p, b, r = u

    J_11 = -λp
    J_12 = -(2*b*kb^2*σp)/(b^2 + kb^2)^2 
    J_13 = -(2*r^3*σp)/(r^2 + kr^2)^2 + (2*r*σp)/(r^2 + kr^2) 

    J_21 = -(2*p*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)^2*(r^2 + kr^2))
    
    bcr_signal_numerator = bcr_max_signal*
        exp(-(t - bcr_max_signal_centered_on_timestep)^2/
        (2*bcr_max_signal_timestep_std^2))  
    bcr_signal_denominator = (sqrt(2*π) * bcr_max_signal_timestep_std)
    bcr_signal_term = bcr_signal_numerator/bcr_signal_denominator 
    dissociation_prod_term = (2*b*kb^2*kp^2*kr^2*σb)/
        ((b^2+kb^2)^2*(p^2+kp^2)*(r^2+kr^2))
    J_22 = -bcr_signal_term - λb - dissociation_prod_term 

    J_23 = -(2*r*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)*(r^2 + kr^2)^2)

    J_31 = 0
    J_32 = 0
    J_33 = -λr - (2*r^3*σr)/((r^2 + kr^2)^2) + (2*r*σr)/(r^2 + kr^2)

    J = [J_11 J_12 J_13
         J_21 J_22 J_23
         J_31 J_32 J_33]

    SMatrix{3, 3}(J)
end 
