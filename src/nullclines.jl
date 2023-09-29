# Functions/models for computing nullclines of equations defined in Martinez2012
# to model germinal cell center regulation

using PlasmaCellDiff
using NaNMath; nm = NaNMath
using UnPack

"""
    blimp1_nullcline(bcl6_level, params::GerminalCenterODEParams)  

```math
\\begin{equation}
\\dot{p} = 0 \\Leftrightarrow n_p(b) = p = \\frac{\\mu_p + \\sigma_p \\frac{k_b^2}{k_b^2 + b^2}}{\\lambda_p}
\\end{equation}
```

Return nullcline for BLIMP1 (`p`).

# References
[1] : Equation S6 from Martinez2012

[2] : https://www.normalesup.org/~doulcier/teaching/modeling/bistable_systems.html
"""
function blimp1_nullcline(bcl6_level, params::GerminalCenterODEParams)
    # parameters 
    @unpack μp, σp, λp, kb = params

    # transcription factor state variables 
    b = bcl6_level

    kb_scaled = dissociation_scaler(kb, b)
    p_nullcline = (μp + σp*kb_scaled)/λp
 
    return p_nullcline
end 

"""
    bcl6_nullcline(bcl6_level, params::GerminalCenterODEParams) 

```math
\\begin{equation}
\\dot{x} = 0 \\Leftrightarrow n_b(b) = \\sqrt{\\frac{k_p^2}{[(\\lambda_b + bcr_0 \\frac{k_b^2}{k_b^2 + b^2})b - \\mu_b]\\frac{k_b^2 + b^2}{\\sigma_b k_b^2}} - k_p^2}
\\end{equation}
```

Return nullcline for BCL6 (`b`). 

NOTE: After discussion with Dr. A. van Kampen, numerically computing
the nullcilnes is also an option (instead of using closed form solution like
the one provided in the current function. The package options are 
nullclines {phaseR} (in CRAN) or homerolled grind.R (university utrecht)

# References
[1] : Equation S7 from Martinez2012

[2] : https://www.normalesup.org/~doulcier/teaching/modeling/bistable_systems.html
"""
function bcl6_nullcline(bcl6_level, params::GerminalCenterODEParams)
    # parameters 
    @unpack kp, kb, σb, μb, λb, bcr0 = params

    b = bcl6_level
  
    # hand solved using bcr0 
    denom = (b*(λb + bcr0*(kb^2/(kb^2 + b^2))) - μb)*(kb^2 + b^2)/(σb*kb^2)
    b_nullcline = nm.sqrt((kp^2/denom) - kp^2)
     
    return b_nullcline
end 
