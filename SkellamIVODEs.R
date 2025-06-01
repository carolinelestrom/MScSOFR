

### E_{tau_7-}[exp(u * X_{tau_8-})]
beta_q_tau8 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-min(T - S, T - jumps[8]) * (1 + x))) - 1) + 
    mu_d * (exp(-nu_D * (-min(T - S, T - jumps[8]) * (1 + x))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau8 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-min(T - S, T - jumps[8]) * (1 + x))) - 1) +
    mu_d * (exp(-nu_D * (-min(T - S, T - jumps[8]) * (1 + x))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_6-}[
###            E_{tau_7-}[exp(u * X_{tau_8-})]]
beta_q_tau7 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[7]) + min(T - S, T - jumps[7]) * (1 + x)))) - 1) +
    mu_d * (exp(-nu_D * (-((S - jumps[7]) + min(T - S, T - jumps[7]) * (1 + x)))) - 1) +
    beta_q_tau8(x, tau = jumps_[8] - jumps_[7])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau7 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[7]) + min(T - S, T - jumps[7]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[7]) + min(T - S, T - jumps[7]) * (1 + x)))) - 1) +
    beta_q_tau8(x, tau = jumps_[8] - jumps_[7])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_5-}[
###            E_{tau_6-}[
###                       E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_q_tau6 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[6]) + min(T - S, T - jumps[6]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[6]) + min(T - S, T - jumps[6]) * (1 + x)))) - 1) + 
    beta_q_tau7(x, tau = jumps_[7] - jumps_[6])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau6 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[6]) + min(T - S, T - jumps[6]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[6]) + min(T - S, T - jumps[6]) * (1 + x)))) - 1) + 
    beta_q_tau7(x, tau = jumps_[7] - jumps_[6])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### E_{tau_4-}[
###            E_{tau_5-}[
###                       E_{tau_6-}[
###                                  E_{tau_7-}[exp(u * X_{tau_8-})]]]
beta_q_tau5 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[5]) + min(T - S, T - jumps[5]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[5]) + min(T - S, T - jumps[5]) * (1 + x)))) - 1) + 
    beta_q_tau6(x, tau = jumps_[6] - jumps_[5])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau5 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[5]) + min(T - S, T - jumps[5]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[5]) + min(T - S, T - jumps[5]) * (1 + x)))) - 1) + 
    beta_q_tau6(x, tau = jumps_[6] - jumps_[5])
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
beta_q_tau4 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[4]) + min(T - S, T - jumps[4]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[4]) + min(T - S, T - jumps[4]) * (1 + x)))) - 1) + 
    beta_q_tau5(x, tau = jumps_[5] - jumps_[4])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau4 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[4]) + min(T - S, T - jumps[4]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[4]) + min(T - S, T - jumps[4]) * (1 + x)))) - 1) + 
    beta_q_tau5(x, tau = jumps_[5] - jumps_[4])
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
beta_q_tau3 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[3]) + min(T - S, T - jumps[3]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[3]) + min(T - S, T - jumps[3]) * (1 + x)))) - 1) + 
    beta_q_tau4(x, tau = jumps_[4] - jumps_[3])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau3 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[3]) + min(T - S, T - jumps[3]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[3]) + min(T - S, T - jumps[3]) * (1 + x)))) - 1) + 
    beta_q_tau4(x, tau = jumps_[4] - jumps_[3])
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
beta_q_tau2 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[2]) + min(T - S, T - jumps[2]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[2]) + min(T - S, T - jumps[2]) * (1 + x)))) - 1) + 
    beta_q_tau3(x, tau = jumps_[3] - jumps_[2])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau2 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[2]) + min(T - S, T - jumps[2]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[2]) + min(T - S, T - jumps[2]) * (1 + x)))) - 1) + 
    beta_q_tau3(x, tau = jumps_[3] - jumps_[2])
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
beta_q_tau1 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[1]) + min(T - S, T - jumps[1]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[1]) + min(T - S, T - jumps[1]) * (1 + x)))) - 1) + 
    beta_q_tau2(x, tau = jumps_[2] - jumps_[1])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_q_tau1 <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-((S - jumps[1]) + min(T - S, T - jumps[1]) * (1 + x)))) - 1) + 
    mu_d * (exp(-nu_D * (-((S - jumps[1]) + min(T - S, T - jumps[1]) * (1 + x)))) - 1) + 
    beta_q_tau2(x, tau = jumps_[2] - jumps_[1])
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}




