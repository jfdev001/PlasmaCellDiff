using Distributions: Normal, pdf
import ConstructionBase: setproperties, setproperties_object
using ConstructionBase
using UnPack


"""
    AbstractGeneRegulationType

Supertype of all concrete implementations used to determine which BCR/CD40
regulation method is called via multiple dispatch.
"""
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
    μx::Float64 = NaN  

    σp::Float64 = 9.0   # Maximum induced transcription rate
    σb::Float64 = 100.0 
    σr::Float64 = 2.6
    σx::Float64 = NaN

    λp::Float64 = 1.0   # Degradation rate 
    λb::Float64 = 1.0 
    λr::Float64 = 1.0
    λx::Float64 = NaN

    kp::Float64 = 1.0   # dissociation constants
    kb::Float64 = 1.0 
    kr::Float64 = 1.0 
    kx::Float64 = NaN

    # Reciprocal function BCR and CD40 regulation parameters
    bcr0::Float64 = 0.05 # Range of BCR-induced degradation of BCL6 in [0, 10]
    cd0::Float64 = 0.015 # Range of CD40-induced transcription of IRF4 in [0, 1] 

    # Gaussian bcr0 and cd0 regulation parameters
    bcr0_max_signal::Float64 = 1
    bcr0_max_signal_centered_on_timestep::Float64 = 45
    bcr0_max_signal_timestep_std::Float64 = 0.1

    cd0_max_signal::Float64 = 1
    cd0_max_signal_centered_on_timestep::Float64 = 60
    cd0_max_signal_timestep_std::Float64 = 0.1

    # For fixed point/bifurcation calculation, the bcr0/cd0 regulatory 
    # mechanism require explicit knowledge about the timestep
    regulation_timestep_t::Float64 = NaN

    # false to eliminate IRF4-mediateded BLIMP1 activation as in fig s3,
    # defaults to true and thus equation S1 is unchanged
    irf4_mediated_blimp1_activation::Bool = true
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
    println(io, ")")
end

"""
    setproperties_object(obj::GerminalCenterODEParams, patch)

Return `GerminalCenterODEParams` constructed with new parameters as specified
by `patch`. Compliant with interface needed for `BifurcationKit.continuation`.

# References:
[1] : [Get parameters of a parametric type](https://stackoverflow.com/questions/35759794/get-parameters-of-a-parametric-type)
"""
function setproperties_object(obj::GerminalCenterODEParams, patch)
    ConstructionBase.check_properties_are_fields(obj)                                               
    nt = getproperties(obj)                                                        
    nt_new = merge(nt, patch)                                                      
    ConstructionBase.check_patch_properties_exist(nt_new, nt, obj, patch)                           
    bcr0_cd0_types = collect(typeof(obj).parameters)
    return GerminalCenterODEParams{bcr0_cd0_types...}(; nt_new...)    
end  

function setproperties(obj::GerminalCenterODEParams, patch::NamedTuple)
    setproperties_object(obj, patch)
end

function Base.getindex(A::GerminalCenterODEParams, inds...) 
    prop_names = propertynames(A)
    prop_values = [getproperty(A, prop_name) for prop_name in prop_names]
    return prop_values[inds...]
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
