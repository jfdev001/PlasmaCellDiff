using UnPack

"""
     check_irf4_bistability_conditions(; μr, cd40, σr, λr, kr)

Prints whether the current IRF4 parameters meet the conditions for bistability.

# References
[1] : Equations S9 - S11 from Martinez2012
"""
function check_irf4_bistability_conditions(; μr, cd40, σr, λr, kr)
    β = (μr + cd40 + σr)/(λr*kr)
    
    if β > sqrt(3)
        println("β > √3" * " => $β > $√3")
    end 
    
    β_polynomial = β^3 - (β^2 - 3)^(3/2) + 9*β 
    β_polynomial_str = "β^3 - (β^2 - 3)^(3/2) + 9*β"
    scaled_transcription_rate_threshold = ((27*σr)/(2*λr*kr))
    scaled_transcription_rate_threshold_str = "((27*σr)/(2*λr*kr))" 

    if β_polynomial < scaled_transcription_rate_threshold 
       println(β_polynomial_str * " < " * 
            scaled_transcription_rate_threshold_str *
            " => " * "$β_polynomial < $scaled_transcription_rate_threshold")
    end 

    if β_polynomial > scaled_transcription_rate_threshold  
       println(β_polynomial_str * " > " * 
            scaled_transcription_rate_threshold_str *
            " => " * "$β_polynomial > $scaled_transcription_rate_threshold")
    end 

    return nothing
end 

function check_irf4_bistability_conditions(params)
    @unpack μr, cd0, σr, λr, kr = params
    check_irf4_bistability_conditions(; μr, cd40 = cd0, σr, λr, kr)
end 
