module PlasmaCellDiff

# parameters
export GerminalCenterODEParams
export Constant, Gaussian

# models
export germinal_center_exit_pathway, germinal_center_exit_pathway_rule
export bcr_subnetwork_rule, bcr_subnetwork
export cd40_subnetwork_rule, cd40_subnetwork
export germinal_center_exit_pathway_jacobian

# gene regulatory signals
export BCR, CD40

# nullclines
export blimp1_nullcline, bcl6_nullcline

# IRF4 bistability functions for figure 2
export check_irf4_bistability_conditions
export compute_irf4_bistability_β

# misc 
export NORMAL_DISTRIBUTION_PEAK_Y_STD_5, NORMAL_DISTRIBUTION_PEAK_Y_STD_1
export bifurcation

include("common.jl")
include("martinez_germinal_center_model.jl")
include("nullclines.jl")
include("jacobian.jl")
include("irf4_bistability.jl")
include("bifurcation.jl")
include("fixedpoints.jl")

end # module PlasmaCellDiff
