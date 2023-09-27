# Kinetic model of the gene regulatory module controlling the germinal 
# center (gc) exit pathway by using ODEs. 
# Implementation based on 
# Martínez MR, et. al. Quantitative modeling of the terminal differentiation of 
# B cells and mechanisms of lymphomagenesis. Proc Natl Acad Sci U S A. 2012 Feb 
# 14; 109(7):2672-7. (aka Martinez2012)

using DynamicalSystems: SVector, SMatrix, CoupledODEs
using Distributions: Normal, pdf
using UnPack

const NORMAL_DISTRIBUTION_PEAK_Y_STD_5::Float64 = 0.08
const NORMAL_DISTRIBUTION_PEAK_Y_STD_1::Float64 = 0.40

"""
    mutable struct GerminalCenterODEParams{BCRMethod, CD40Method}

Parameters used for dynamical systems simulating germinal center regulation.
`BCRMethod` and `CD40Method` are symbols in 
`[:constant, :gaussian, :reciprocal]` that determine how CD40 and BCR
are modelled. Note that `:reciprocal` uses equations S4 and S5 defined 
explicitly in Martinez2012.
"""
@kwdef mutable struct GerminalCenterODEParams{BCRMethod, CD40Method}
    μp::Float64 = 10e-6 # Basal transcription rate
    μb::Float64 = 2.0 
    μr::Float64 = 0.1  

    σp::Float64 = 9.0   # Maximum induced transcription rate
    σb::Float64 = 100.0 
    σr::Float64 = 2.6

    kp::Float64 = 1.0 
    kb::Float64 = 1.0 
    kr::Float64 = 1.0 

    λp::Float64 = 1.0   # Degradation rate 
    λb::Float64 = 1.0 
    λr::Float64 = 1.0

    # BCR and CD40 gene regulation as constant in time
    bcr_constant::Float64 = 15.0 # from figure S1
    cd40_constant::Float64 = NaN

    # Reciprocal function BCR and CD40 regulation parameters
    bcr0::Float64 = 0.05 # Range of BCR-induced degradation of BCL6 in [0, 10]
    cd0::Float64 = 0.015 # Range of CD40-induced transcription of IRF4 in [0, 1] 
    C0::Float64 = 10e-8  # Appears to be unused

    # Gaussian BCR and CD40 regulation parameters
    # determined experimentally in 
    # `notebooks/plot_martinez_germinal_center_gaussian_trajectory.ipynb`
    bcr_max_signal::Float64 = 0.1875
    bcr_max_signal_centered_on_timestep::Float64 = 50
    bcr_max_signal_timestep_std::Float64 = 1

    cd40_max_signal::Float64 = 0.0375
    cd40_max_signal_centered_on_timestep::Float64 = 60
    cd40_max_signal_timestep_std::Float64 = 1
end

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
    dissociation_scaler(k, ui)

Scales dissociation constant k associated with transcription factor uᵢ.
"""
dissociation_scaler(k, ui) = k^2 / (k^2 + ui^2)

"""
    transcription_factor_scaler(k, ui)

Scales protein level for transcription factor ui using dissociation constant
k associated with uᵢ.
"""
transcription_factor_scaler(k, ui) = ui^2 / (k^2 + ui^2)

"""
    BCR(; bcr0, kb, b)

Return BCR gene regulatory signal.

# References
[1] : Equation S4 from Martinez2012
"""
BCR(; bcr0, kb, b) = bcr0*dissociation_scaler(kb, b)

function BCR(u, params::GerminalCenterODEParams{:constant, T}, t) where T
    @unpack bcr_constant = params
    return bcr_constant 
end 

function BCR(u, params::GerminalCenterODEParams{:gaussian, T}, t) where T 
    # gaussian regulation parameters 
    @unpack bcr_max_signal = params
    @unpack bcr_max_signal_centered_on_timestep = params
    @unpack bcr_max_signal_timestep_std = params

    return gaussian_regulatory_signal(;
        peak = bcr_max_signal, 
        μ = bcr_max_signal_centered_on_timestep,
        σ = bcr_max_signal_timestep_std,
        t = t)
end 

function BCR(u, params::GerminalCenterODEParams{:reciprocal, T}, t) where T
    @unpack bcr0, kb = params
    p, b, r = u
    return BCR(; bcr0, kb, b)
end 

"""
    CD40(; cd0, kb, b)

Return CD40 gene regulatory signal.

# References
[1] : Equation S5 from Martinez2012
"""
CD40(; cd0, kb, b) = cd0*dissociation_scaler(kb, b)

function CD40(u, params::GerminalCenterODEParams{T, :constant}, t) where T
    @unpack cd40_constant = params
    return cd40_constant
end

function CD40(u, params::GerminalCenterODEParams{T, :gaussian}, t) where T
    @unpack cd40_max_signal = params
    @unpack cd40_max_signal_centered_on_timestep = params
    @unpack cd40_max_signal_timestep_std = params
 
    return gaussian_regulatory_signal(; 
        peak = cd40_max_signal,
        μ = cd40_max_signal_centered_on_timestep,
        σ = cd40_max_signal_timestep_std,
        t = t)
end 

function CD40(u, params::GerminalCenterODEParams{T, :reciprocal}, t) where T
    @unpack cd0, kb = params
    p, b, r = u
    return CD40(; cd0, kb, b)
end 

"""
    gaussian_regulatory_signal(; peak, μ, σ, t)

Return BCR/CD40 regulatory signal from evaluating the PDF of the normal 
distribution with the desired parameters `μ` and `σ` at point `t` scaled by 
`peak`.

NOTE: Repeatedly instantiating `Normal` is likely not efficient, but the 
alternative is passing it as a parameter to the function defining the ODE
system (`germinal_center_gaussian_exit_pathway`) which is known to be
inefficient. Perhaps `Benchmark`ing this would be interesting, though is
clearly not necessary for the scope of this project.
"""
gaussian_regulatory_signal(; peak, μ, σ, t) = peak*pdf(Normal(μ, σ), t)

"""
    irf4_bistability(; μr, cd40, σr, λr, k) 

Return bistability constraint for IRF4 written as a function of relevant IRF4
kinetic parameters.

NOTE: Maybe altering σr is what leads to the Figures 2 (β ∈ {1.5, 1.8, 2})

# References
[1] : Equation S9 (i.e., β = ...) from Martinez2012
"""
irf4_bistability(; μr, cd40, σr, λr, kr) = (μr + cd40 + σr)/(λr*kr)


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
