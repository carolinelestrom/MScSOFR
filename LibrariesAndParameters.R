setwd("~/Documents/KU/MSc/FinalRCode")
####################################################################################################################
####################################################################################################################
#------------------------------------- Libraries -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
library(dplyr) ### Data Manipulation
library(ggplot2) ### Plots
library(deSolve) ### Differential Equations
library(lubridate) ### Date Formats
library(pracma) ### Integration
library(deSolve) ### Ordinary Differential Equations
library(rmutil) ### Ordinary Differential Equations

####################################################################################################################
####################################################################################################################
#------------------------------------ Parameters -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### IORB
r0 <- 0.03

### Continuous SOFR Spread
kappa_s_jump <- 0.0442
theta_s_jump <- 0.1056
sigma_s_jump <- 0.1
lambda_s_jump <- -0.4783
s0 <- 0.003

### Normal Unscheduled Jump Sizes
nu_P <- 0.428
mu_P <- -0.0075
sigma_P <- 0.01

### Mixed-Exponential Unscheduled Jump Sizes
nu_P <- 0.428
pu_J <- res_CDF$pu_J #c()
qd_J <- res_CDF$qd_J #c()
lambda_J <- res_CDF$lambda_J #c()
eta_J <- res_CDF$eta_J #c()
q_J <- res_CDF$q_J #c()
theta_J <- res_CDF$theta_J #c()

### Normal Scheduled Jump Sizes
omega <- 0.01
gamma_Q <- 0
Gamma_Q <- c(1,0)
kappa_r <- 1.6988
kappa_xi <- kappa_r
kappa_theta <- 1.7046
theta_theta <- -0.0003
sigma_xi <- 0.1
sigma_theta <- 0.1
rho <- -0.9447
xi0 <- theta_theta ### -0.0003
theta0 <- theta_theta ### -0.0003

### Skellam Scheduled Jump Sizes
nu_D <- 0.0008
mu_u <- 3333
mu_d <- 3336.75
gamma_Qu <- 0
gamma_Qd <- 0
GammaQu <- mu_u
GammaQd <- mu_d
gamma_Pu <- 
gamma_Pd <- 
Gamma_Pu <- 
Gamma_Pd <- 
lambda_u <-
lambda_d <- 
theta_zeta <- 0.1
kappa_zeta <- 1.6988
sigma_zeta <- 0.1
zeta0 <- theta_zeta

nu_D * (GammaQu * MeanXSkellam(0, 0.1, zeta0) - GammaQd * MeanXSkellam(0, 0.1, zeta0))

nu_D * (GammaQu * theta_zeta - GammaQd * theta_zeta)
