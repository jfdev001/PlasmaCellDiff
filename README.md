# PlasmaCellDiff

Repository for the period 1 2023 masters course Bioinformatics 1 at the University of Amsterdam.

## Progress Meeting 2 (2023-09-22)

Things to include in order of priority

* [X] CDR40 and BCR signals (use trajectory outputs and compute per timestep)
* [X] Read Chaper 8 on Cell for understanding derivations of Martinez2012 model
* [X] Interactive trajectory/3D trajectory 
* Bifurcation analysis
* Sensitivity analysis
* Chaos or other entropy?
* Critical points

## Progress Meeting 3 (2023-10-06)

* Tasks (dynamic):
    * [X] Figure S1 (assigned)
        * [X] Nullcline analysis (Ch. 8 Cell and Ch.2 JuliaNLD)
            * BCR signaling might be increased or decreased by modulating the
            magnitude of the gaussian peak
        * [X] Bifurcation/hysteresis plots
            * BCR ?= bcr0, so vary this and then change $\lambda_b$ and $\sigma_b$,
            degradation and transcription parameters, respectively. 
    * [X] Figure 2 (assigned)
        * [X] Email discussion suggests impossible
        * [X] Show pitchfork '>' and saddle point bifurcation 'S' anyways
        * Generally could not reproduce
        * [] Could attempt naive fixed point computation for different
        initial conditions and only the single moving parameter cd0 ...
    * [X] Figure 3 (assigned)
        * BCR behavior replicated but CD40 not
    * [] Analytical bifurcation of gaussian model and determination
        of why the behavior was not replicated 
    * [] Sensitivity analysis (should be pretty mindless, call global
        sensitivity on dynamical system -- with constant bcr0/cd0
        and could possibly see how this is with gaussian bcr0/cd0
    * [] Figure S3
        * This is just the gaussian model and the thick lines in the figure
        are probably replicable by simply r^2/kr^2 + r^2 term from eq S1
        (i.e., "elimination of IRF4-mediated BLIMP1 activation"

## On Fixed Point Analysis for Figure 3 and Possibly Figure S5

$$
\frac{dp}{dt} = \mu_p + \sigma_p \frac{k_b^2}{k_b^2 + b^2} + \sigma_p \frac{r^2}{k_r^2 + r^2} - \lambda_p p
$$

$$
\frac{db}{dt} = \mu_b + \sigma_b \frac{k_p^2}{k_p^2 + p^2}\frac{k_b^2}{k_b^2 + b^2}\frac{k_r^2}{k_r^2 + r^2} - (\lambda_b + BCR)b
$$

$$
\frac{dr}{dt} = \mu_r + \sigma_r \frac{r^2}{k_r^2 + r^2} + CD40 - \lambda_r r 
$$

$$
BCR = BCR_{max} \frac{1}{\sigma_{BCR}\sqrt{2\pi}}\exp[-\frac{1}{2}(\frac{t - \mu_{BCR}}{\sigma_{BCR}})^2]\frac{k_b^2}{k_b^2 + b^2}
$$

$$
CD40 = CD40_{max}\frac{1}{\sigma_{CD40}\sqrt{2\pi}}\exp[-\frac{1}{2}(\frac{t - \mu_{CD40}}{\sigma_{CD40}})^2]\frac{k_b^2}{k_b^2 + b^2}
$$

# References

[Good Example of  Scientific Project in Julia: HighDimensionalComplexityEntropy](https://github.com/ikottlarz/HighDimensionalComplexityEntropy)

Datseris G. and Parlitz U. Chapter 1 - 4 from Nonlinear Dynamics: A Concise 
Introduction Interlaced with Code (2022). 
[github](https://github.com/JuliaDynamics/NonlinearDynamicsTextbook/tree/master). 
note: null clines, bifurcation analysis, stability analysis, lyapunov.

Numerical bifurcation diagrams and bistable systems. Lecture notes from
computational biology at École normale supérieure de Paris. 
url: http://www.normalesup.org/~doulcier/teaching/modeling/ (2018)
