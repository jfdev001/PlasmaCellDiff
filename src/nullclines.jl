# Functions/models for computing nullclines of equations defined in Martinez2012
# to model germinal cell center regulation

using UnPack
include("common.jl")

"""
    blimp1_nullcline(u, params::GerminalCenterODEParams)  

```math
\\begin{equation}
\\dot{p} = 0 \\Leftrightarrow  n_p(b) = \\sigma_p + \\sigma_p \\frac{k_b^2}{k_b^2 + b^2}
\\end{equation}
```

Return nullcline for BLIMP1 (`p`).

# References
[1] : Equation S6 from Martinez2012

[2] : https://www.normalesup.org/~doulcier/teaching/modeling/bistable_systems.html
"""
function blimp1_nullcline(u, params::GerminalCenterODEParams)     
    # parameters 
    @unpack μp, μb = params
    @unpack σp, σb = params
    @unpack kp, kb = params
    @unpack λp, λb = params

    # transcription factor state variables 
    p, b = u

    # pdot = 0, therefore p = Integral[pdot] = 0, and n_p(b) = ...
    kb_scaled = dissociation_scaler(kb, b)
    p_nullcline = μp + σp*kb_scaled
 
    return p_nullcline
end 

function bcl6_nullcline(u, params::GerminalCenterODEParams) 
    # parameters 
    @unpack μp, μb = params
    @unpack σp, σb = params
    @unpack kp, kb = params
    @unpack λp, λb = params

    # transcription factor state variables 
    p, b = u
    
    
end 
