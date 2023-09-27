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
    μₚ::Float64 = 10e-6, # Basal transcription rate
    μb::Float64 = 2.0, 
    μᵣ::Float64 = 0.1,  

    σₚ::Float64 = 9.0,   # Maximum induced transcription rate
    σb::Float64 = 100.0, 
    σᵣ::Float64 = 2.6,

    kₚ::Float64 = 1.0, 
    kb::Float64 = 1.0, 
    kᵣ::Float64 = 1.0, 

    λₚ::Float64 = 1.0,   # Degradation rate 
    λb::Float64 = 1.0, 
    λᵣ::Float64 = 1.0,

    # BCR and CD40 gene regulation as constant in time
    bcr_constant::Float64 = 15.0, # from figure S1
    cd40_constant::Float64 = NaN,

    # Reciprocal function BCR and CD40 regulation parameters
    bcr₀::Float64 = 0.05, # Range of BCR-induced degradation of BCL6 in [0, 10]
    cd₀::Float64 = 0.015, # Range of CD40-induced transcription of IRF4 in [0, 1] 
    C₀::Float64 = 10e-8,  # Appears to be unused

    # Gaussian BCR and CD40 regulation parameters
    # determined experimentally in 
    # `notebooks/plot_martinez_germinal_center_gaussian_trajectory.ipynb`
    bcr_max_signal::Float64 = 0.1875,
    bcr_max_signal_centered_on_timestep::Float64 = 50, 
    bcr_max_signal_timestep_std::Float64 = 1,

    cd40_max_signal::Float64 = 0.0375,
    cd40_max_signal_centered_on_timestep::Float64 = 60,
    cd40_max_signal_timestep_std::Float64 = 1
end

"""
    germinal_center_exit_pathway_rule(u, params::GerminalCenterODEParams, t)

Return rule for gene regulatory module controlling germinal center exit pathway
with coupled BCR and CD40 regulatory signals.

# Arguments 
- `u`: State vector of dynamical system. The elements of this vector represent
    protein levels of `p` (BLIMP1), `b` (BCL6), and `r` (IRF4).
- `p::GerminalCenterODEParams`: Typed mutable struct whose parametrized
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
    @unpack μₚ, μb, μᵣ = params
    @unpack σₚ, σb, σᵣ = params
    @unpack kₚ, kb, kᵣ = params
    @unpack λₚ, λb, λᵣ = params

    # transcription factor state variables 
    p, b, r = u

    # compute scaled dissociation constants and protein levels
    kp_scaled = dissociation_scaler(kₚ, p)
    kb_scaled = dissociation_scaler(kb, b)
    kr_scaled = dissociation_scaler(kᵣ, r)
    r_scaled = transcription_factor_scaler(kᵣ, r)

    # regulatory signals
    bcr = BCR(; u, p = params, t)
    cd40 = CD40(; u, p = params, t)

    # system describing evolution of transcription factors in germinal center
    pdot = μₚ + σₚ*kb_scaled + σₚ*r_scaled - λₚ*p
    bdot = μb + σb*kp_scaled*kb_scaled*kr_scaled - (λb + bcr)*b
    rdot = μᵣ + σᵣ*r_scaled + cd40 - λᵣ*r

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
    dissociation_scaler(k, uᵢ)

Scales dissociation constant k associated with transcription factor uᵢ.
"""
dissociation_scaler(k, uᵢ) = k^2 / (k^2 + uᵢ^2)

"""
    transcription_factor_scaler(k, uᵢ)

Scales protein level for transcription factor uᵢ using dissociation constant
k associated with uᵢ.
"""
transcription_factor_scaler(k, uᵢ) = uᵢ^2 / (k^2 + uᵢ^2)

"""
    BCR(; bcr₀, kb, b)

Return BCR gene regulatory signal.

# References
[1] : Equation S4 from Martinez2012
"""
BCR(; bcr₀, kb, b) = bcr₀*dissociation_scaler(kb, b)

function BCR(; u, p::GerminalCenterODEParams{:constant, T}, t) where T
    @unpack bcr_constant = p
    return bcr_constant 
end 

function BCR(; u, p::GerminalCenterODEParams{:gaussian, T}, t) where T 
    # gaussian regulation parameters 
    @unpack bcr_max_signal = p
    @unpack bcr_max_signal_centered_on_timestep = p
    @unpack bcr_max_signal_timestep_std = p

    return gaussian_regulatory_signal(;
        peak = bcr_max_signal, 
        μ = bcr_max_signal_centered_on_timestep,
        σ = bcr_max_signal_timestep_std,
        t = t)
end 

function BCR(; u, p::GerminalCenterODEParams{:reciprocal, T}, t) where T
    @unpack bcr₀, kb
    p, b, r = u
    return BCR(; bcr₀, kb, b)
end 

"""
    CD40(; cd₀, kb, b)

Return CD40 gene regulatory signal.

# References
[1] : Equation S5 from Martinez2012
"""
CD40(; cd₀, kb, b) = cd₀*dissociation_scaler(kb, b)

function CD40(u, p::GerminalCenterODEParams{T, :constant}, t) where T
    @unpack cd40_constant = p
    return cd40_constant
end

function CD40(u, p::GerminalCenterODEParams{T, :gaussian}, t) where T
    @unpack cd40_max_signal = p
    @unpack cd40_max_signal_centered_on_timestep = p
    @unpack cd40_max_signal_timestep_std = p
 
    return gaussian_regulatory_signal(; 
        peak = cd40_max_signal,
        μ = cd40_max_signal_centered_on_timestep,
        σ = cd40_max_signal_timestep_std,
        t = t)
end 

function CD40(u, p::GerminalCenterODEParams{T, :reciprocal}, t) where T
    @unpack cd₀, kb
    p, b, r = u
    return CD40(; cd₀, kb, b)
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
    irf4_bistability(; μᵣ, cd40, σᵣ, λᵣ, k) 

Return bistability constraint for IRF4 written as a function of relevant IRF4
kinetic parameters.

NOTE: Maybe altering σᵣ is what leads to the Figures 2 (β ∈ {1.5, 1.8, 2})

# References
[1] : Equation S9 (i.e., β = ...) from Martinez2012
"""
irf4_bistability(; μᵣ, cd40, σᵣ, λᵣ, kᵣ) = (μᵣ + cd40 + σᵣ)/(λᵣ*kᵣ)


"""
    germinal_center_exit_pathway_jacobian(
        u, p::GerminalCenterODEParams{:reciprocal, :reciprocal}, t)

Return Jacobian of `germinal_center_exit_pathway_rule` w/ reciprocal BCR/CD40.

See `notebooks/germinal_cell_jacobian.pdf` for the source that generated the
symbolic Jacobian.

TODO: Check if both are requried to be reciprocal for this to be true.
"""
function germinal_center_exit_pathway_jacobian(
    u, p::GerminalCenterODEParams{:reciprocal, :reciprocal}, t)
    @unpack μₚ, μb, μᵣ = p
    @unpack σₚ, σb, σᵣ = p
    @unpack kₚ, kb, kᵣ = p
    @unpack λₚ, λb, λᵣ = p
    @unpack bcr₀, cd₀, C₀ = p

    p, b, r = u

    J_11 = -λₚ
    J_12 = -(2*b*kb^2*σₚ)/(b^2 + kb^2)^2 
    J_13 = -(2*r^3*σₚ)/(r^2 + kᵣ^2)^2 + (2*r*σₚ)/(r^2 + kᵣ^2) 

    J_21 = -(2*p*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)^2*(r^2 + kᵣ^2))
    J_22 = (2*b^2*bcr₀*kb^2)/(b^2+kb^2)^2 - (bcr₀*kb^2)/(b^2 + kb^2) - 
        λb - (2*b*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)^2*(p^2 + kₚ^2)*(r^2 + kᵣ^2))
    J_23 = -(2*r*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)*(r^2 + kᵣ^2)^2)
 
    J_31 = 0
    J_32 = -(2*b*cd₀*kb^2)/(b^2 + kb^2)^2
    J_33 = -λᵣ - (2*r^3*σᵣ)/(r^2 + kᵣ^2)^2 + (2*r*σᵣ)/(r^2 + kᵣ^2)

    J = [J_11 J_12 J_13
         J_21 J_22 J_23
         J_31 J_32 J_33]

    return SMatrix{3, 3}(J)
end 

"""
    germinal_center_exit_pathway_jacobian(
        u, p::GerminalCenterODEParams{:gaussian, :gaussian}, t) 

Return Jacobian of `germinal_center_exit_pathway_rule` w/ gaussian BCR/CD40.

See `notebooks/germinal_cell_gaussian_jacobian.pdf` for the source that 
generated the symbolic Jacobian.

TODO: check if both are required to be gaussian for this Jacobian to be correct.
"""
function germinal_center_exit_pathway_jacobian(
    u, p::GerminalCenterODEParams{:gaussian, :gaussian}, t) 
    # parameters 
    @unpack μₚ, μb, μᵣ = p
    @unpack σₚ, σb, σᵣ = p
    @unpack kₚ, kb, kᵣ = p
    @unpack λₚ, λb, λᵣ = p
    @unpack bcr₀, cd₀, C₀ = p
   
    # gaussian regulation parameters 
    @unpack bcr_max_signal = p
    @unpack bcr_max_signal_centered_on_timestep = p
    @unpack bcr_max_signal_timestep_std = p

    @unpack cd40_max_signal = p
    @unpack cd40_max_signal_centered_on_timestep = p
    @unpack cd40_max_signal_timestep_std = p

    p, b, r = u

    J_11 = -λₚ
    J_12 = -(2*b*kb^2*σₚ)/(b^2 + kb^2)^2 
    J_13 = -(2*r^3*σₚ)/(r^2 + kᵣ^2)^2 + (2*r*σₚ)/(r^2 + kᵣ^2) 

    J_21 = -(2*p*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)^2*(r^2 + kᵣ^2))
    
    bcr_signal_numerator = bcr_max_signal*
        exp(-(t - bcr_max_signal_centered_on_timestep)^2/
        (2*bcr_max_signal_timestep_std^2))  
    bcr_signal_denominator = (sqrt(2*π) * bcr_max_signal_timestep_std)
    bcr_signal_term = bcr_signal_numerator/bcr_signal_denominator 
    dissociation_prod_term = (2*b*kb^2*kₚ^2*kᵣ^2*σb)/
        ((b^2+kb^2)^2*(p^2+kₚ^2)*(r^2+kᵣ^2))
    J_22 = -bcr_signal_term - λb - dissociation_prod_term 

    J_23 = -(2*r*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)*(r^2 + kᵣ^2)^2)

    J_31 = 0
    J_32 = 0
    J_33 = -λᵣ - (2*r^3*σᵣ)/((r^2 + kᵣ^2)^2) + (2*r*σᵣ)/(r^2 + kᵣ^2)

    J = [J_11 J_12 J_13
         J_21 J_22 J_23
         J_31 J_32 J_33]

    SMatrix{3, 3}(J)
end 
