# PlasmaCellDiff

Codebase implementing the germinal center exit pathway kinetic model, emphasizing bifurcation analysis through simulation, as introduced by Martinez et. al. (2012). Part of the masters course Bioinformatics 1 during period 1, 2023 at the University of Amsterdam.

## Gaussian Signaled Germinal Cell Center Model

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

M. R. Martínez et al., Proceedings of the National Academy of Sciences, vol. 109, no. 7, pp. 2672–2677, 2012. doi:10.1073/pnas.1113019109

[Good Example of  Scientific Project in Julia: HighDimensionalComplexityEntropy](https://github.com/ikottlarz/HighDimensionalComplexityEntropy)

Datseris G. and Parlitz U. Chapter 1 - 4 from "Nonlinear Dynamics: A Concise Introduction Interlaced with Code"(2022). [github](https://github.com/JuliaDynamics/NonlinearDynamicsTextbook/tree/master). note: null clines, bifurcation analysis, stability analysis, lyapunov exponents.

Numerical bifurcation diagrams and bistable systems. Lecture notes from
computational biology at École normale supérieure de Paris.
url: http://www.normalesup.org/~doulcier/teaching/modeling/ (2018)
