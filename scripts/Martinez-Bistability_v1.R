###################################################
#
# Martinez model of plasma cell differentation
# 
# Reference
# Mart?nez et al (2012) Quantitative modeling of the terminal 
# differentiation of B cells and mechanisms of lymphomagenesis. 
# PNAS, 109:2672-7.
#
# Bistability analysis using the R package phaseR
#
# Bistability analyses of the BCR submodel:
#   Because we want to understand the dynamics of
#   BCR signaling when decoupled from CD40 signaling, we assume
#   that the levels of IRF4 are much smaller than kr and, therefore,
#   the protein levels do not change substantially during the early
#   response phase.
#   Bistability condition: sigma_b / lambda_b > 20
#
# Calculation of nullclines
#
# Implemented by Antoine van Kampen
# Bioinformatics Laboratory
# Amsterdam University Medical Centers (Location AMC)
#
###################################################

library(deSolve)
library(phaseR)

#
#clear workspace
#
rm(list=ls())

# # Getting the path of your current open file
# # This prevents to manually specify the directory with setwd()
# current_path = rstudioapi::getActiveDocumentContext()$path 
# setwd(dirname(current_path ))
# getwd()

parameters = c(
    mu_p      = 10^-6,  # Basal transcription rate  BLIMP
    mu_b      = 2,      # Basal transcription rate BCL6
    mu_r      = 0.1,    # Basal transcription rate IRF4
    sigma_p   = 9,      # Maximum induced transcription rate BLIMP
    sigma_b   = 100,    # Maximum induced transcription rate BLC6
    sigma_r   = 2.6,    # Maximum induced transcription rate IRF4
    k_p       = 1,      # Dissociation constant BLIMP
    k_b       = 1,      # Dissociation constant BLC6
    k_r       = 1,      # Dissociation constant IRF4
    lambda_p  = 1,      # Degration rate BLIMP
    lambda_b  = 1,      # Degration rate BLC6
    lambda_r  = 1)      # Degration rate IRF4


############################
# Derivatives
############################
Martinez = function(t, y, parms) {

  with(as.list(c(y,parms)), {
  
  a_bcr=10
  a_cd=5
  s=0.1
  mean_bcr=45
  mean_cd=60
  bcr0=a_bcr*exp(-(t-mean_bcr)^2/s^2)  #gaussian distributed bcr0 (signal spike)
  cd0=a_cd*exp(-(t-mean_cd)^2/s^2)     #gaussian distributed bcr0 (signal spike)

  BCR = bcr0*(k_b^2/(k_b^2+BCL6^2))
  CD40 = cd0*(k_b^2/(k_b^2+BCL6^2))

  BCR_t[i,]<<-c(t,BCR)
  CD40_t[i,]<<-c(t,CD40)
  i<<-i+1

  dBLIMP1 = mu_p + sigma_p*(k_b^2/(k_b^2+BCL6^2)) + sigma_p*(IRF4^2/(k_r^2+IRF4^2)) - lambda_p*BLIMP1
  dBCL6 =   mu_b + sigma_b*(k_p^2/(k_p^2+BLIMP1^2))*(k_b^2/(k_b^2+BCL6^2))*(k_r^2/(k_r^2+IRF4^2)) - (lambda_b+BCR)*BCL6
  dIRF4 =   mu_r + sigma_r*(IRF4^2/(k_r^2+IRF4^2)) + CD40 - lambda_r*IRF4

  l=list(c(dBLIMP1,dBCL6,dIRF4))
  return(l) })
}

#
# First plot the solution curves
#
set.seed(34)
i <<- 1
BCR_t<<-data.frame(matrix(ncol=2,dimnames=list(c(),c("Time","BCR"))))
CD40_t<<-data.frame(matrix(ncol=2,dimnames=list(c(),c("Time","CD50"))))

# initial conditions for each clone
yini = c(BLIMP1 = 0, BCL6 = 5, IRF4 = 0)
times =  seq(from = 0, to = 200, by = 0.5) # time span in days; 
out = ode(y = yini, times = times, func = Martinez, parms = parameters, method = "lsoda")

matplot(out[,1],out[,-1],main="Transcription factor levels",xlab='Time',ylab="Level",lwd=2,lty=1,type='l')
points(BCR_t[,1],5*BCR_t[,2],type='l',col='orange') 
points(CD40_t[,1],5*CD40_t[,2],type='l',col='blue')
legend("topleft",legend=c("BCR","CD40"),col=c("orange","blue"),lty=c(1,1))


##
# Next, make phase diagram of the BCR network with phaseR
##

BCR_Network=function (t, y, parameters) {
  bcr0     = 11 #see Figure S1 Martinez; I used 11 instead of 15
  BCL6     = y[1]
  BLIMP1   = y[2]

  mu_b     = parameters['mu_b']
  mu_p     = parameters['mu_p']
  sigma_b  = parameters['sigma_b']
  sigma_p  = parameters['sigma_p']
  k_p      = parameters['k_p']
  k_b      = parameters['k_b']
  lambda_b = parameters['lambda_b']
  lambda_p = parameters['lambda_p']
  
  dy     = numeric(2)

  dy[1] = mu_b + sigma_b*(k_p^2/(k_p^2+BLIMP1^2))*(k_b^2/(k_b^2+BCL6^2)) - 
          (lambda_b+bcr0*(k_b^2/(k_b^2+BCL6^2)))*BCL6
  dy[2] = mu_p + sigma_p*(k_b^2/(k_b^2+BCL6^2)) - lambda_p*BLIMP1 
  list(dy)
}

BCR_Network.flowField = flowField(BCR_Network, 
                                  xlim = c(0,6), ylim = c(0,10), 
                                  parameters=parameters, 
                                  points = 20, col="grey50",
                                  system="two.dim", 
                                  xlab="BCL6",ylab="BLIMP1",
                                  add = FALSE)

y0 = matrix(c(0, 1, 0, 4, 1, 0, 2, 6, 4, 4), ncol = 2, nrow = 5, byrow = TRUE)

BCR_Network.nullclines <- nullclines(BCR_Network, 
                                     xlim = c(0,6), ylim = c(0,10),col=c('blue','red'),
                                     parameters = parameters, 
                                     points = 100,
                                     add.legend = FALSE)

legend("topright",legend=c("BCL6","BLIMP1"),
       text.col=c("red","blue"),col=c("red","blue"),cex=0.9,lty=1)

BCR_Network.trajectory <- trajectory(BCR_Network, y0 = y0, tlim = c(0,10),
                                     parameters = parameters)

#Evaluate stability (Jacobian) for points in the phase plane
stability(BCR_Network, parameters=parameters, 
          system="two.dim",
          ystar = c(0.5,8.2))

stability(BCR_Network, parameters=parameters, 
          system="two.dim",
          ystar = c(1.9,2))

