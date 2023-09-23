# Kinetic model of the gene regulatory module controlling the germinal 
# center (gc) exit pathway by using ODEs. 
# Implementation based on 
# Martínez MR, et. al. Quantitative modeling of the terminal differentiation of 
# B cells and mechanisms of lymphomagenesis. Proc Natl Acad Sci U S A. 2012 Feb 
# 14; 109(7):2672-7. (aka Martinez2012)

using DynamicalSystems: SVector, SMatrix, CoupledODEs
using Distributions: Normal, pdf
using UnPack

const NORMAL_DISTRIBUTION_PEAK_Y = 0.08
germinal_center_ode_params = Dict{Symbol, Float64}(
    :μₚ => 10e-6, # Basal transcription rate
    :μb => 2.0, 
    :μᵣ => 0.1,  

    :σₚ => 9.0,   # Maximum induced transcription rate
    :σb => 100.0, 
    :σᵣ => 2.6,

    :kₚ => 1.0, 
    :kb => 1.0, 
    :kᵣ => 1.0, 

    :λₚ => 1.0,   # Degradation rate 
    :λb => 1.0, 
    :λᵣ => 1.0,

    :bcr₀ => 0.05, # Range of BCR-induced degradation of BCL6 in [0, 10]
    :cd₀ => 0.015,  # Range of CD40-induced transcription of IRF4 in [0, 1] 
    :C₀ => 10e-8, 

    # gaussian BCR and CD40 regulation parameters
    # determined experimentally in 
    # `notebooks/plot_martinez_germinal_center_gaussian_trajectory.ipynb`
    :bcr_max_signal => 0.015 / NORMAL_DISTRIBUTION_PEAK_Y,
    :bcr_max_signal_centered_on_timestep => 50, 
    :bcr_max_signal_timestep_std => 1,

    :cd40_max_signal => 0.003 / NORMAL_DISTRIBUTION_PEAK_Y,
    :cd40_max_signal_centered_on_timestep => 60,
    :cd40_max_signal_timestep_std => 1,
)

"""
    germinal_center_exit_pathway_rule(u, p, t)

Return rule for gene regulatory module controlling germinal center exit pathway
with coupled BCR and CD40 regulatory signals.

# Arguments 
- `u`: State vector of dynamical system. The elements of this vector represent
    protein levels of `p` (BLIMP1), `b` (BCL6), and `r` (IRF4).
- `p`: Parameters for the model.
- `t`: Current time for numerical integration. No need to pass an argument here
    since numerical integration is handled with builtin solvers.

# References
[1] : Equations S1 - S3 from Martinez2012
"""
function germinal_center_exit_pathway_rule(u, p, t)
    # parameters 
    @unpack μₚ, μb, μᵣ = p
    @unpack σₚ, σb, σᵣ = p
    @unpack kₚ, kb, kᵣ = p
    @unpack λₚ, λb, λᵣ = p
    @unpack bcr₀, cd₀, C₀ = p

    # transcription factor state variables 
    p, b, r = u

    # compute scaled dissociation constants and protein levels
    kp_scaled = dissociation_scaler(kₚ, p)
    kb_scaled = dissociation_scaler(kb, b)
    kr_scaled = dissociation_scaler(kᵣ, r)
    r_scaled = transcription_factor_scaler(kᵣ, r)

    # regulatory signals
    bcr = BCR(; bcr₀, kb, b)
    cd40 = CD40(; cd₀, kb, b)

    # system describing evolution of transcription factors in germinal center
    pdot = μₚ + σₚ*kb_scaled + σₚ*r_scaled - λₚ*p
    bdot = μb + σb*kp_scaled*kb_scaled*kr_scaled - (λb + bcr)*b
    rdot = μᵣ + σᵣ*r_scaled + cd40 - λᵣ*r

    return SVector(pdot, bdot, rdot)
end

"""
    germinal_center_gaussian_exit_pathway_rule(u, p, t)

Return rule for gene regulatory module controlling germinal center exit pathway
with coupled gaussian BCR and CD40 gene regulatory signals.

Implements CD40 and BCR gene regulation as gaussian distribution centered
on a different timesteps as well as having different maximum signals
"""
function germinal_center_gaussian_exit_pathway_rule(u, p, t)
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

    # transcription factor state variables 
    p, b, r = u

    # compute scaled dissociation constants and protein levels
    kp_scaled = dissociation_scaler(kₚ, p)
    kb_scaled = dissociation_scaler(kb, b)
    kr_scaled = dissociation_scaler(kᵣ, r)
    r_scaled = transcription_factor_scaler(kᵣ, r)

    # regulatory signals
    bcr = gaussian_regulatory_signal(; 
        peak = bcr_max_signal, 
        μ = bcr_max_signal_centered_on_timestep, 
        σ = bcr_max_signal_timestep_std,
        t = t)
    cd40 = gaussian_regulatory_signal(; 
        peak = cd40_max_signal,
        μ = cd40_max_signal_centered_on_timestep,
        σ = cd40_max_signal_timestep_std,
        t = t)

    # system describing evolution of transcription factors in germinal center
    pdot = μₚ + σₚ*kb_scaled + σₚ*r_scaled - λₚ*p
    bdot = μb + σb*kp_scaled*kb_scaled*kr_scaled - (λb + bcr)*b
    rdot = μᵣ + σᵣ*r_scaled + cd40 - λᵣ*r

    return SVector(pdot, bdot, rdot)
end 

"""
    germinal_center_exit_pathway(u0, params) 

Return CoupledODEs model for gene regulatory module controlling germinal center 
exit pathway with coupled BCR and CD40 regulatory signals. 
"""
function germinal_center_exit_pathway(u0, params) 
    return CoupledODEs(germinal_center_exit_pathway_rule, u0, params)
end 

"""
    germinal_center_gaussian_exit_pathway(u0, params) 

Return CoupledODEs model for gene regulatory module controlling germinal center 
exit pathway with gaussian BCR and CD40 regulatory signals. 

TODO: Should a seed be set?? If so where? Does it make sense to have
the calling of the random distribution func in the ds rule??
"""
function germinal_center_gaussian_exit_pathway(u0, params) 
    return CoupledODEs(germinal_center_gaussian_exit_pathway_rule, u0, params)
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

"""
    CD40(; cd₀, kb, b)

Return CD40 gene regulatory signal.

# References
[1] : Equation S5 from Martinez2012
"""
CD40(; cd₀, kb, b) = cd₀*dissociation_scaler(kb, b)

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
    germinal_center_exit_pathway_jacobian(u, p, t)

Return Jacobian of `germinal_center_exit_pathway_rule`.

See `notebooks/germinal_cell_jacobian.pdf` for the source that generated the
symbolic Jacobian.
"""
function germinal_center_exit_pathway_jacobian(u, p, t)
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
    germinal_center_gaussian_exit_pathway_jacobian(u, p, t)

Return Jacobian of `germinal_center_gaussian_exit_pathway_rule`.

See `notebooks/germinal_cell_gaussian_jacobian.pdf` for the source that 
generated the symbolic Jacobian.
"""
function germinal_center_gaussian_exit_pathway_jacobian(u, p, t)
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


