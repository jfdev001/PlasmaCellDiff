"""
     bifurcation(f, bifurcation_param, total_time, transient_time)

Compute bifurcations of the ODE system `f` by varying the `bifurcation_param`
from `p_min` to `p_max` and evolving the system for some `transient_time` 
before saving (equlibrium) states for a system evolving for some `total_time`. 

# Return
the bifurcation parameters and equilibrium states of `f`.
"""
function bifurcation(
    f, bifurcation_param, p_min, p_max, total_time, transient_time)

end 
