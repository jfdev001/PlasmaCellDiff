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

TODO: This could be made more efficient by using design decision
inspired by Chaitanya Kumar (and his suggestions about solvers used
in DifferentialEquations.jl) by employing compile time dispatch as well as
using a non mutable struct (as suggested by 
https://stackoverflow.com/questions/77188675/overloading-methods-for-parametrized-struct
```julia
abstract type AbstractParamType end

struct Exponential <: AbstractParamType end
struct Gaussian <: AbstractParamType end
struct Constant <: AbstractParamType end
struct Reciprocal <: AbstractParamType end

@kwdef mutable struct Params{P1<:AbstractParamType, P2<:AbstractParamType}

compute_bcr(;u, params::Params{Constant, <:AbstractParamType}, t) = 15
```
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

function Base.show(io::IO, z::GerminalCenterODEParams)
    propnames = propertynames(z)
    propvalues = [getproperty(z, propname_i) for propname_i in propnames]
    println("$(typeof(z))(")
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

