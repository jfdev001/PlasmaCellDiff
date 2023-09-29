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

* [] Figure S1
    * [X] Nullcline analysis (Ch. 8 Cell and Ch.2 JuliaNLD)
        * BCR signaling might be increased or decreased by modulating the
        magnitude of the gaussian peak
    * [] Bifurcation/hysteresis plots
        * BCR ?= bcr0, so vary this and then change $\lambda_b$ and $\sigma_b$,
        degradation and transcription parameters, respectively. 
* [] Figure 2
* [] FIgure 3

# References

[Good Example of  Scientific Project in Julia: HighDimensionalComplexityEntropy](https://github.com/ikottlarz/HighDimensionalComplexityEntropy)

Datseris G. and Parlitz U. Chapter 1 - 4 from Nonlinear Dynamics: A Concise 
Introduction Interlaced with Code (2022). 
[github](https://github.com/JuliaDynamics/NonlinearDynamicsTextbook/tree/master). 
note: null clines, bifurcation analysis, stability analysis, lyapunov.

Numerical bifurcation diagrams and bistable systems. Lecture notes from
computational biology at École normale supérieure de Paris. 
url: http://www.normalesup.org/~doulcier/teaching/modeling/ (2018)
