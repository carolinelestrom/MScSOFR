

### E_{tau_7-}[exp(u * X_{tau_8-})]
beta_tau8 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[8]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[8]))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau8 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[8]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[8]))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_6-}[
###            E_{tau_7-}[exp(u * X_{tau_8-})]]
beta_tau7 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[7]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[7]))) - 1) +
    beta_tau8(tau = jumps_[8] - jumps_[7])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau7 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[7]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[7]))) - 1) +
    beta_tau8(tau = jumps_[8] - jumps_[7])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_5-}[
###            E_{tau_6-}[
###                       E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau6 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[6]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[6]))) - 1) +
    beta_tau7(tau = jumps_[7] - jumps_[6])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau6 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[6]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[6]))) - 1) +
    beta_tau7(tau = jumps_[7] - jumps_[6])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_4-}[
###            E_{tau_5-}[
###                       E_{tau_6-}[
###                                  E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau5 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[5]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[5]))) - 1) +
    beta_tau6(tau = jumps_[6] - jumps_[5])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau5 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[5]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[5]))) - 1) +
    beta_tau6(tau = jumps_[6] - jumps_[5])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_3-}[
###            E_{tau_4-}[
###                       E_{tau_5-}[
###                                  E_{tau_6-}[
###                                             E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau4 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[4]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[4]))) - 1) +
    beta_tau5(tau = jumps_[5] - jumps_[4])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau4 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[4]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[4]))) - 1) +
    beta_tau5(tau = jumps_[5] - jumps_[4])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_2-}[
###            E_{tau_3-}[
###                       E_{tau_4-}[
###                                  E_{tau_5-}[
###                                             E_{tau_6-}[
###                                                        E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau3 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[3]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[3]))) - 1) +
    beta_tau4(tau = jumps_[4] - jumps_[3])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau3 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[3]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[3]))) - 1) +
    beta_tau4(tau = jumps_[4] - jumps_[3])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_1-}[
###            E_{tau_2-}[
###                       E_{tau_3-}[
###                                  E_{tau_4-}[
###                                             E_{tau_5-}[
###                                                        E_{tau_6-}[
###                                                                   E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau2 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[2]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[2]))) - 1) +
    beta_tau3(tau = jumps_[3] - jumps_[2])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau2 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[2]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[2]))) - 1) +
    beta_tau3(tau = jumps_[3] - jumps_[2])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_t[
###     E_{tau_1-}[
###                E_{tau_2-}[
###                           E_{tau_3-}[
###                                      E_{tau_4-}[
###                                                 E_{tau_5-}[
###                                                            E_{tau_6-}[
###                                                                       E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_tau1 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[1]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[1]))) - 1) +
    beta_tau2(tau = jumps_[2] - jumps_[1])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_tau1 <- function(tau){
  u <- mu_u * (exp(nu_D * (-(T - jumps[1]))) - 1) + mu_d * (exp(-nu_D * (-(T - jumps[1]))) - 1) +
    beta_tau2(tau = jumps_[2] - jumps_[1])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}



exp(alpha_tau8(tau = jumps_[8] - jumps_[7]) + 
      alpha_tau7(tau = jumps_[7] - jumps_[6]) + 
      alpha_tau6(tau = jumps_[6] - jumps_[5]) + 
      alpha_tau5(tau = jumps_[5] - jumps_[4]) + 
      alpha_tau4(tau = jumps_[4] - jumps_[3]) + 
      alpha_tau3(tau = jumps_[3] - jumps_[2]) + 
      alpha_tau2(tau = jumps_[2] - jumps_[1]) + 
      alpha_tau1(tau = jumps_[1] - 0) + 
      beta_tau1(tau = jumps_[1] - 0) * zeta0)





