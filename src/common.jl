using Distributions: Normal, pdf
using UnPack

abstract type AbstractGeneRegulationType end

struct Gaussian <: AbstractGeneRegulationType end
struct Constant <: AbstractGeneRegulationType end

"""
    struct GerminalCenterODEParams{
        Bcr0Method <: AbstractGeneRegulationType, 
        Cd0Method <: AbstractGeneRegulationType} 

Parameters used for dynamical systems simulating germinal center regulation.
`Bcr0Method` and `Cd0Method` denote whether the gene regulation signals
propagated by BCR and CD40 will be modified by either a constant bcr0/cd0
or a gaussian bcr0/cd0 (see [`AbstractGeneRegulationType`](@ref) for valid 
types).

NOTE: Could make the gaussian distribution for bcr0 and cd0 fields of this
struct; however, mixed parameter type might cause performance issues
for dynamical systems... current strategy is to just instantiate 
a `Normal()` distribution each time see [`gaussian_regulatory_signal`](@ref)

# References
[1] : Table S1 from Martinez2012

[2]: [Stack Overflow: On abstract parameter type design](https://stackoverflow.com/questions/77188675/overloading-methods-for-parametrized-struct)
"""
@kwdef struct GerminalCenterODEParams{
    Bcr0Method <: AbstractGeneRegulationType, 
    Cd0Method <: AbstractGeneRegulationType} 
    # parameters
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

    # Reciprocal function BCR and CD40 regulation parameters
    bcr0::Float64 = 0.05 # Range of BCR-induced degradation of BCL6 in [0, 10]
    cd0::Float64 = 0.015 # Range of CD40-induced transcription of IRF4 in [0, 1] 
    C0::Float64 = 10e-8  # Appears to be unused

    # Gaussian bcr0 and cd0 regulation parameters
    # determined experimentally in 
    bcr0_max_signal::Float64 = 1
    bcr0_max_signal_centered_on_timestep::Float64 = 45
    bcr0_max_signal_timestep_std::Float64 = 0.1
    # TODO: 
    # bcr_gaussian = Normal(...) # save peak for actual function call

    cd0_max_signal::Float64 = 1
    cd0_max_signal_centered_on_timestep::Float64 = 60
    cd0_max_signal_timestep_std::Float64 = 0.1
    # cd40_gaussian = Normal(...)
end

function Base.show(io::IO, z::GerminalCenterODEParams)
    propnames = propertynames(z)
    propvalues = [getproperty(z, propname_i) for propname_i in propnames]
    println(io, "$(typeof(z))(")
    for i in 1:length(propnames)
        pname = propnames[i]
        pvalue = propvalues[i]
        println(io, "  $pname = $pvalue")
    end 
    println(")")
end 

"""
    dissociation_scaler(k, ui)

Scales dissociation constant k associated with transcription factor uᵢ.
"""
dissociation_scaler(k, ui) = k^2 / (k^2 + ui^2)

"""
    transcription_factor_scaler(k, ui)

Scales protein level for transcription factor ui using dissociation constant
k associated with ui.
"""
transcription_factor_scaler(k, ui) = ui^2 / (k^2 + ui^2)

"""
    BCR(; bcr0, kb, b)

Return BCR gene regulatory signal.

# References
[1] : Equation S4 from Martinez2012
"""
BCR(; bcr0, kb, b) = bcr0*dissociation_scaler(kb, b)

function BCR(
    u, 
    params::GerminalCenterODEParams{Gaussian, <:AbstractGeneRegulationType}, 
    t)  
    # gaussian regulation parameters 
    @unpack bcr0_max_signal = params
    @unpack bcr0_max_signal_centered_on_timestep = params
    @unpack bcr0_max_signal_timestep_std = params
    @unpack kb = params

    b = u[2]

    gaussian_bcr0 = gaussian_regulatory_signal(;
        peak = bcr0_max_signal, 
        μ = bcr0_max_signal_centered_on_timestep,
        σ = bcr0_max_signal_timestep_std,
        t = t)

    bcr = BCR(; bcr0 = gaussian_bcr0, kb, b)

    return bcr
end 

function BCR(
    u, 
    params::GerminalCenterODEParams{Constant, <:AbstractGeneRegulationType}, 
    t) 
    @unpack bcr0, kb = params
    b = u[2]
    return BCR(; bcr0, kb, b)
end 

"""
    CD40(; cd0, kb, b)

Return CD40 gene regulatory signal.

# References
[1] : Equation S5 from Martinez2012
"""
CD40(; cd0, kb, b) = cd0*dissociation_scaler(kb, b)

function CD40(
    u, 
    params::GerminalCenterODEParams{<:AbstractGeneRegulationType, Gaussian}, 
    t)
    @unpack cd0_max_signal = params
    @unpack cd0_max_signal_centered_on_timestep = params
    @unpack cd0_max_signal_timestep_std = params
    @unpack kb = params
    
    b = u[2]

    gaussian_cd0 = gaussian_regulatory_signal(; 
        peak = cd0_max_signal,
        μ = cd0_max_signal_centered_on_timestep,
        σ = cd0_max_signal_timestep_std,
        t = t)

    cd40 = CD40(; cd0 = gaussian_cd0, kb, b)

    return cd40
end 

function CD40(
    u, 
    params::GerminalCenterODEParams{<:AbstractGeneRegulationType, Constant}, 
    t)
    @unpack cd0, kb = params
    b = u[2]
    return CD40(; cd0, kb, b)
end 

"""
    gaussian_regulatory_signal(; peak, μ, σ, t)

Return BCR/CD40 regulatory signal from evaluating the PDF of the normal 
distribution with the desired parameters `μ` and `σ` at point `t` scaled by 
`peak`.

NOTE: Repeatedly instantiating `Normal` is likely not efficient, but the 
alternative is passing it as a parameter to the function defining the ODE
system.
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

