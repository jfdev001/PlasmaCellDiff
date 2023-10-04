using PlasmaCellDiff
using DynamicalSystems: SMatrix
using UnPack

"""
    germinal_center_exit_pathway_jacobian(
        u, params::GerminalCenterODEParams{Constant, Constant}, t)

Return Jacobian of `germinal_center_exit_pathway_rule` w/ constant cd0/bcr0.

See `notebooks/germinal_cell_jacobian.pdf` for the source that generated the
symbolic Jacobian.
"""
function germinal_center_exit_pathway_jacobian(
    u, params::GerminalCenterODEParams{Constant, Constant}, t)
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
        u, params::GerminalCenterODEParams{Gaussian, Gaussian}, t) 

Return Jacobian of `germinal_center_exit_pathway_rule` w/ gaussian BCR/CD40.
"""
function germinal_center_exit_pathway_jacobian(
    u, params::GerminalCenterODEParams{Gaussian, Gaussian}, t) 

    # parameters 
    @unpack μp, μb, μr = params 
    @unpack σp, σb, σr = params 
    @unpack kp, kb, kr = params 
    @unpack λp, λb, λr = params 
    @unpack bcr0, cd0, C0 = params
   
    # gaussian regulation parameters 
    @unpack bcr0_max_signal = params
    @unpack bcr0_max_signal_centered_on_timestep = params
    @unpack bcr0_max_signal_timestep_std = params
    bcr0m = bcr0_max_signal 
    bcr0μ = bcr0_max_signal_centered_on_timestep 
    bcr0σ = bcr0_max_signal_timestep_std
    
    @unpack cd0_max_signal = params
    @unpack cd0_max_signal_centered_on_timestep = params
    @unpack cd0_max_signal_timestep_std = params    
    cd0m = cd0_max_signal 
    cd0μ = cd0_max_signal_centered_on_timestep 
    cd0σ = cd0_max_signal_timestep_std

    # state variables
    p, b, r = u

    # pdot/dui
    J_11 = -λp
    J_12 = -(2*b*kb^2*σp)/(b^2 + kb^2)^2 
    J_13 = -(2*r^3*σp)/(r^2 + kr^2)^2 + (2*r*σp)/(r^2 + kr^2) 

    # bdot/dui
    J_21 = -(2*p*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)^2*(r^2 + kr^2)) 
    
    bcr0_exp = exp(-(-bcr0μ + t)^2/(2*bcr0σ^2))
    bcr_term1 = (sqrt(2)*b^2*bcr0m*bcr0_exp*kb^2)/(bcr0σ*(b^2 + kb^2)^2*sqrt(π))
    bcr_term2 = (bcr0m*bcr0_exp*kb^2)/(bcr0σ*(b^2 + kb^2)*sqrt(2*π))
    J_22 = bcr_term1 - bcr_term2 - λb - (2*kb*kb^2*kp^2*kr^2*σb)/
        ((b^2 + kb^2)^2*(kp^2 + p^2)*(kr^2 + r^2))

    J_23 = -(2*r*kb^2*kp^2*kr^2*σb)/((b^2 + kb^2)*(p^2 + kp^2)*(r^2 + kr^2)^2)

    # rdot/dui
    J_31 = 0
    
    cd0_exp = exp(-(-cd0μ + t)^2/(2*cd0σ^2))
    J_32 = - (sqrt(2)*b*cd0m*cd0_exp*kb^2)/(cd0σ*(b^2 + kb^2)^2*sqrt(π)) 

    J_33 = -λr - (2*r^3*σr)/((r^2 + kr^2)^2) + (2*r*σr)/(r^2 + kr^2)

    J = [J_11 J_12 J_13
         J_21 J_22 J_23
         J_31 J_32 J_33]

    SMatrix{3, 3}(J)
end 
