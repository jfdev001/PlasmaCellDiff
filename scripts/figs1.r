## Plotting nullclines (fig s1 martinez2012) using phaseR

# https://hluebbering.github.io/phase-planes/
library(knitr)
library(phaseR)
library(deSolve)
library(graphics)
library(latex2exp)

# Figure options
knitr::opts_chunk$set(
  echo = FALSE, out.width = 400, fig.align = "center", message = FALSE)

# Open png device
# https://stackoverflow.com/questions/7144118/how-to-save-a-plot-as-image-on-the-disk
png("figures/figs1_nullclines_phaseR.png")

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

bcr_subnetwork_flowField <-flowField(bcr_subnetwork,
     xlim = c(0.0, 5.5),
     ylim = c(0., 5.5),
     add  = FALSE,
     main = "{phaseR} Computed Nullclines for BCR Subnetwork",
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