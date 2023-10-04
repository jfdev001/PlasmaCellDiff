using DynamicalSystems: trajectory

""" 
    function bifurcation(
        ds::Function, u0, params,
        bifurcation_param::Symbol, p_min, p_max, length_param_space,
        total_time, transient_time) 

NOTE: just remove this... isn't functionally different from
ChaosTools.orbitdiagram

Compute bifurcations of the dynamical system `ds` by varying the 
`bifurcation_param` from `p_min` to `p_max` with a total parameter space 
`length_param_space` then evolve the system from the initial conditions `u0` 
for a `total_time` and keeping only the states **after** the `transient_time`.

# Arguments
- `ds::Function`: Function that returns a dynamical system of the form 
        `ds(u0, params)`.

# Return
The bifurcation parameters and equilibrium states of `f`.

# References
[1] : https://discourse.julialang.org/t/construct-namedtuple-dynamically/15394/7

[2] : https://scicomp.stackexchange.com/questions/41994/python-bifurcation-diagram-of-seasonally-forced-epidemiological-models
"""
function bifurcation(
    ds::Function, u0, params,
    bifurcation_param::Symbol, p_min, p_max, length_param_space,
    total_time, transient_time) 

    # vectors to store repeated param for plotting as well as the steady states
    bifurcation_plot_params = []
    bifurcation_plot_trajectories = []

    # iterate through param space, evolve system, and save bifurcation plot data
    param_space = LinRange(p_min, p_max, length_param_space)
    for param_value in param_space
        bifurcation_param_nt = (bifurcation_param, )
        param_value_nt = (param_value,)
        new_param = (;zip(bifurcation_param_nt, param_value_nt)...)
        params = setproperties(params, new_param)
        myds = ds(u0, params)
        X, t = trajectory(myds, total_time; Ttr = transient_time)
        repeated_param_value = repeat([param_value], length(X))
        push!(bifurcation_plot_params, repeated_param_value)
        push!(bifurcation_plot_trajectories, X)
    end

    return bifurcation_plot_params, bifurcation_plot_trajectories
end 
