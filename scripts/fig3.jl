using PlasmaCellDiff
using DynamicalSystems
using Plots
using UnPack

# figure 3A with abrogated CD40 signaling
u0 = [0.2, 5.0, 0.2]
bcr0_max_signal = 155
params = GerminalCenterODEParams{Gaussian, Gaussian}(
    regulation_timestep_t = 20,
    
    cd0_max_signal = 0.0, # 0.80 * bcr0_max_signal
    cd0_max_signal_centered_on_timestep = 50,
    cd0_max_signal_timestep_std = 5,
    
    bcr0_max_signal = bcr0_max_signal,
    bcr0_max_signal_centered_on_timestep = 35,
    bcr0_max_signal_timestep_std = 5)

ds = germinal_center_exit_pathway(u0, params)

# compute fixed points
p_domain = b_domain = r_domain = interval(0, 10)
state_domain = p_domain × b_domain × r_domain
fp, eigs, stable = fixedpoints(
    ds, state_domain, germinal_center_exit_pathway_jacobian)
