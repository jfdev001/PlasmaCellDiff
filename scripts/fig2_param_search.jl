using PlasmaCellDiff
using UnPack
using BifurcationKit
using Plots
using Dates
using Formatting 
Plots.default(show=false) # prevent showing plot on every call

raw"""
    plot_fig2(
        ; folder = \"tmpfigs\",
        u0 = [0.0], 
        cd0s=[0.0], 
        μrs=0:0.01:0.3,  
        σrs=1.6:0.01:2.8, 
        λrs=0.80:0.01:1.6, 
        krs=0.80:0.01:1.6)

Plots the desired parameter space for figure 2 CD40 subnetwork. 

TODO: Should add a sampling parameter to downsample the parameter
space.
"""
function plot_fig2(
    ; folder = "tmpfigs",
    u0 = [0.0], 
    cd0s=[0.0], 
    μrs=0:0.01:0.2,  
    σrs=1.6:0.01:2.8, 
    λrs=0.80:0.01:1.6, 
    krs=0.80:0.01:1.6,
    report_every_n_iters=1000,)

    @assert isdir(folder) "folder must exist"

    # datetime and open error file
    dtime = Dates.format(now(), "yyyymmdd_H-M-S")

    # Log number of iterations needed
    niters = length(cd0s)*length(μrs)*length(σrs)*length(λrs)*length(krs)
    println("niters = " * format(niters, commas = true))
    println("Sleeping for 5 seconds to give time to cancel this")
    sleep(5)
    println("Beginning iteration")

    # iterate through parameters and plot bifurcationd diagrams
    cnter = 1
    for cd0 in cd0s,  λr in λrs, kr in krs, σr in σrs, μr in μrs,
        params = GerminalCenterODEParams{Constant, Constant}(;
            cd0,
            μr,
            σr,
            λr,
            kr)

        irf4_β = compute_irf4_bistability_β(; μr, cd40 = cd0, σr, λr, kr)

        prob = BifurcationProblem(
            cd40_subnetwork_rule,
            u0,
            params,
            (@lens _.cd0),
            record_from_solution = (u, params) ->(r = u[1]))

        # parameters for continuation
        opts_br = ContinuationPar(
            p_min = 0.0, 
            p_max = 0.3,
            max_steps = 10_000,
            ds = 0.02,
            dsmax = 0.03,)

        
        # initialize name for plot to save
        params_str = "cd0-$(cd0)_mur-$(μr)"*
            "_sigmar-$(σr)_lambdar-$(λr)_kr-$(kr)"

        plot_fname = "$(dtime)_fig2_$(params_str)"

        # set a flag that checks whether these parameters have already been
        # explored
        already_plotted = false

        # continuation of equilibria
        try
            if !already_plotted 
                br = BifurcationKit.continuation(
                    prob, 
                    MoorePenrose(tangent = PALC()),
                    bothside = true,
                    opts_br; normC = norminf)

                bifurcplot = BifurcationKit.plot(
                    br,
                    ylims = (0, maximum(br.x)),
                    legend = false,
                    title = "(β = $(irf4_β),\n" *
                    " cd0 = $(cd0)," *
                    " μr = $(μr)," *
                    " σr = $(σr)" *
                    " λr = $(λr)," *
                    " kr = $(kr))",
                    titlefont = 10)

                fpath = joinpath(folder, plot_fname)
                savefig(fpath * ".png")            
            end 
        catch e
            # pass on convergence issues
        end

        cnter += 1

        if cnter % report_every_n_iters == 0
            println("PROGRESS: ", "$(cnter)/$(niters)")
        end
    end 
end 
