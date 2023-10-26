## Plotting nullclines (fig s1 martinez2012) using phaseR
# NOTE: This is basically just a copy paste of an old tutorial for which
# one of the packages `captioner` may or may not be in CRAN anymore.
# I didn't want to waste time fiddling with other plot libraries... so sorry
# for the overhead of installing `captioner`

# https://hluebbering.github.io/phase-planes/
library(knitr)
library(phaseR)
library(deSolve)
library(graphics)
# https://cran.r-project.org/web/packages/captioner/index.html
# https://www.dataquest.io/blog/install-package-r/ 
library(captioner)
library(latex2exp)

# Figure options
knitr::opts_chunk$set(
  echo = FALSE, out.width = 400, fig.align = "center", message = FALSE)
fig_nums <- captioner(prefix = "Figure")

# Martinez2012 nullclines 
dissociation_scaler <- function(k, ui){
  (k^2)/(k^2 + ui^2)
}

bcr_subnetwork <- function(t, y, parameters){
    # BLIMP1 parameters
    kp <- 1
    lambda_p <- 1
    sigma_p <- 5
    mu_p <- 10e-6
    
    # BCL6 parameters
    lambda_b <- 1
    bcr0 <- 15
    mu_b <- 1
    sigma_b <- 100
    kb <- 1
  
    # transcription factor state variables
    b <- y[1]
    p <- y[2]
    
    # ODE system
    kb_scaled = dissociation_scaler(kb, b)
    kp_scaled = dissociation_scaler(kp, p)
    bcr <- bcr0*kb_scaled
    
    pdot <- mu_p + sigma_p*kb_scaled - lambda_p*p
    bdot <- mu_b + sigma_b*kp_scaled*kb_scaled - (lambda_b + bcr)*b
    
    du <- c(pdot, bdot)
    list(du)
}

body_cap7 <- fig_nums(
  name = "phase7", 
  caption = "Phase Plane Analysis. Determine the steady-state solution by the 
  nullclines' intersection for the Martinez2012 BCR Subnetwork.")
bcr_subnetwork_flowField <-flowField(bcr_subnetwork,
                                     xlim = c(0.0, 5.5),
                                     ylim = c(0., 5.5),
                                     add  = FALSE,
                                     ylab = "BLIMP1",
                                     xlab= "BCL6",
                                     frac=1)
grid()
bcr_subnetwork_nullclines <- nullclines(
  bcr_subnetwork,
  xlim = c(0.0, 5.5),
  ylim = c(0.0, 5.5),
  lty = 2, lwd = 2,
  col=c("blue","red"))
