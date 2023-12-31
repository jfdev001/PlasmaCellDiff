# `ChaosTools.fixedpoints` uses `t=0.0` during root finding. I defined
# functions here that override the existing root finding helper functions
# so that a non-zero timestep invokes gaussian signaling of BCR/CD40
using PlasmaCellDiff
import ChaosTools: to_root_f, to_root_J, dynamic_rule, fixedpoints
using ChaosTools
using LinearAlgebra
using IntervalRootFinding

# Regulation timestep will also affect dynamical system, not just jacobian
function to_root_f(
    ds::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing)
    return u -> dynamic_rule(ds)(u, p, p.regulation_timestep_t)
end 

function to_root_J(
    Jf, ::CoupledODEs, p::GerminalCenterODEParams{Gaussian, Gaussian}, ::Nothing) 
    return u -> Jf(u, p, p.regulation_timestep_t)
end

"""
    fixedpoints(
            ds::DynamicalSystem, box, regulation_timestep_property::Symbol,
            J = nothing;
            method = IntervalRootFinding.Krawczyk, tol = 1e-15, warn = true,
            order = nothing)

Return fixed points, eigenvalues, and stability of dynamical system
that has some time dependent forcing/regulation function for which the argument
to such a function is in the parameters of `ds` and the symbol corresponding
to the the timestep is `regulation_timestep_property`.
"""
function fixedpoints(
        ds::DynamicalSystem, box, regulation_timestep_property::Symbol,
        J = nothing;
        method = IntervalRootFinding.Krawczyk, tol = 1e-15, warn = true,
        order = nothing)
    if isinplace(ds)
        error("`fixedpoints` currently works only for out-of-place dynamical systems.")
    end
    # Jacobian: copy code from `DynamicalSystemsBase`
    f = dynamic_rule(ds)
    p = current_parameters(ds)
    regulation_timestep_t = getproperty(p, regulation_timestep_property)

    # estimate jacobian is probably best in this case rather
    # than explicit form...
    if isnothing(J)
        # does the regulation timestep here in the fuunction even
        # make a difference??
        Jf(u, p, t) = DynamicalSystemsBase.ForwardDiff.jacobian(
            x -> f(x, p, regulation_timestep_t), u)
    else
        Jf = J
    end

    # Find roots via IntervalRootFinding.jl
    fun = to_root_f(ds, p, order)
    jac = to_root_J(Jf, ds, p, order)
    r = IntervalRootFinding.roots(fun, jac, box, method, tol)
    D = dimension(ds)
    fp = ChaosTools.roots_to_dataset(r, D, warn)

    # Find eigenvalues and stability (taking regulation timestep into account)
    eigs = Vector{Vector{Complex{Float64}}}(undef, length(fp))
    Jm = zeros(dimension(ds), dimension(ds)) # `eigvals` doesn't work with `SMatrix`
    for (i, u) in enumerate(fp)
        # notice that we use the "pure" jacobian, no -u!
        Jm .= Jf(u, p, regulation_timestep_t)
        eigs[i] = LinearAlgebra.eigvals(Array(Jm))
    end
    stable = Bool[ChaosTools.isstable(ds, e) for e in eigs]
    return fp, eigs, stable
end
