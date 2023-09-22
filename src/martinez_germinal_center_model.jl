# Kinetic model of the gene regulatory module controlling the germinal 
# center (gc) exit pathway by using ODEs. 
# Implementation based on 
# Martínez MR, et. al. Quantitative modeling of the terminal differentiation of 
# B cells and mechanisms of lymphomagenesis. Proc Natl Acad Sci U S A. 2012 Feb 
# 14; 109(7):2672-7. (aka Martinez2012)

using DynamicalSystems
using UnPack

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

    :bcr₀ => NaN, # Range of BCR-induced degradation of BCL6 in [0, 10]
    :cd₀ => NaN,  # Range of CD40-induced transcription of IRF4 in [0, 1] 
    :C₀ => 10e-8, 
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
    # TODO: Make stochastic???
    bcr = BCR(; bcr₀, kb, b)
    cd40 = CD40(; cd₀, kb, b)

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
    bcr_subnetwork(u, p, t)

Dynamics of BCR signaling decoupled from CD40 signaling. 

# References
[1] : Equations S3, S6, and S7 from Martinez2012
"""
function bcr_subnetwork_rule(u, p, t)
    @unpack μₚ, μb, μᵣ = p
    @unpack σₚ, σb, σᵣ = p
    @unpack kₚ, kb, kᵣ = p
    @unpack λₚ, λb, λᵣ = p
    @unpack bcr₀, cd₀, C₀ = p

    p, b, r = u 

    # modified state variables for decoupling BCR from CD40 
    pdot = μₚ + σₚ*(kb^2 / (kb^2 + b^2)) - λₚ*p
    bdot = μb + σb*(kₚ^2 / (kₚ^2 + b^2))*(kb^2 / (kb^2 + b^2)) - 
        (λb + bcr₀*(kb^2 / (kb^2 + b^2)))

    # IRF4 state 
    rdot = μᵣ + σᵣ*(r^2 / (kᵣ^2 + r^2)) + cd40 - λᵣ*r

    return SVector(pdot, bdot, rdot)
end 

"""
    cd40_subnetwork(u, p, t)

Same system as germinal_center_exit_pathway???
"""
function cd40_subnetwork_rule(u, p, t)
    @unpack μₚ, μb, μᵣ = p
    @unpack σₚ, σb, σᵣ = p
    @unpack kₚ, kb, kᵣ = p
    @unpack λₚ, λb, λᵣ = p
    @unpack bcr₀, cd₀, C₀ = p

    p, b, r = u 
    cd40 = CD40(; cd₀, kb, b)

    # modified state variable for absence of signals 
    rdot = nothing 

    throw("notimplemented")
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
These should be gaussian peaks at certain time points see figure 3

peak centered on some time point (should be a parameter)
"""
BCR(; bcr₀, kb, b) = bcr₀*dissociation_scaler(kb, b)

CD40(; cd₀, kb, b) = cd₀*dissociation_scaler(kb, b)

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

    J = zeros(3, 3)

    J[1, 1] = -λₚ
    J[1, 2] = -(2*b*kb^2*σₚ)/(b^2 + kb^2)^2 
    J[1, 3] = -(2*r^3*σₚ)/(r^2 + kᵣ^2)^2 + (2*r*σₚ)/(r^2 + kᵣ^2) 

    J[2, 1] = -(2*p*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)^2*(r^2 + kᵣ^2))
    J[2, 2] = (2*b^2*bcr₀*kb^2)/(b^2+kb^2)^2 - (bcr₀*kb^2)/(b^2 + kb^2) - 
        λb - (2*b*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)^2*(p^2 + kₚ^2)*(r^2 + kᵣ^2))
    J[2, 3] = -(2*r*kb^2*kₚ^2*kᵣ^2*σb)/((b^2 + kb^2)*(p^2 + kₚ^2)*(r^2 + kᵣ^2)^2)
 
    J[3, 1] = 0
    J[3, 2] = -(2*b*cd₀*kb^2)/(b^2 + kb^2)^2
    J[3, 3] = -λᵣ - (2*r^3*σᵣ)/(r^2 + kᵣ^2)^2 + (2*r*σᵣ)/(r^2 + kᵣ^2)

    return SMatrix{3, 3}(J)
end 


