####################################################################################################################
####################################################################################################################
#------------------------------------ Settings ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Parameters
set.seed(NULL) #3333333
n_sim <- 10000
T <- 6/12
dt <- 1 / 360 #0.01
n_steps <- T / dt


jump1 <- 28 / 360
jump1_ <- 27 / 360
jump2 <- (2 * 30 + 18) / 360
jump2_ <- (2 * 30 + 17) / 360
jump3 <- (4 * 30 + 6) / 360
jump3_ <- (4 * 30 + 5) / 360
jump4 <- (5 * 30 + 17) / 360
jump4_ <- (5 * 30 + 16) / 360
jump5 <- (6 * 30 + 29) / 360
jump5_ <- (6 * 30 + 28) / 360
jump6 <- (8 * 30 + 16) / 360
jump6_ <- (8 * 30 + 15) / 360
jump7 <- (9 * 30 + 28) / 360
jump7_ <- (9 * 30 + 27) / 360
jump8 <- (11 * 30 + 9) / 360
jump8_ <- (11 * 30 + 8) / 360
jumps <- c(jump1, jump2, jump3, jump4, jump5, jump6, jump7, jump8)
jumps_ <- c(jump1_, jump2_, jump3_, jump4_, jump5_, jump6_, jump7_, jump8_)



####################################################################################################################
####################################################################################################################
#------------------------------------ Functions --------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Stochastic SOFR Spread with No Jumps
SimSOFRSpread <- function() {
  s <- numeric(n_steps + 1)
  s[1] <- s0
  for (i in 2:(n_steps + 1)) {
    #print(i)
    #s[i] <- s[i - 1] + kappa_s_jump * (theta_s_jump - s[i - 1]) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
    s[i] <- s[i - 1] * exp(-kappa_s_jump * dt) +
      theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
      sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
  }
  return(sum(r0 + s) * dt)
}

### Stochastic SOFR Spread + Normal Scheduled Jumps + Normal Unscheduled Jumps
SimNN <- function() {
  s <- numeric(n_steps + 1)
  s[1] <- s0
  r <- numeric(n_steps + 1)
  r[1] <- r0
  for (i in 2:(n_steps + 1)) {
    t <- i * dt
    #print(i)
    #s[i] <- s[i - 1] + kappa_s_jump * (theta_s_jump - s[i - 1]) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
    s[i] <- s[i - 1] * exp(-kappa_s_jump * dt) +
      theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
      sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
    
    ### Unscheduled jumps
    N_P <- rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * 1/360)
    ### Simulate jump sizes Z_1,Z_2,...
    if (N_P > 0) {
      Z_Normal <- rnorm(N_P, mean = mu_P, sd = sigma_P)
    } else {
      Z_Normal <- 0
    }
    ### Update the jump process
    J_P_Normal <- sum(Z_Normal) ### Compute compound Poisson sum
    
    ### Scheduled jumps
    dW_xi <- rnorm(1)
    dW_theta <- rnorm(1)
    #thetam <- thetam + kappa_theta * (theta_theta - thetam) * dt + sqrt(dt) * sigma_theta * (rho * dW_xi + sqrt(1 - rho^2) * dW_theta)
    thetam <- thetam * exp(-kappa_theta * dt) +
      theta_theta * (1 - exp(-kappa_theta * dt)) +
      sigma_theta * sqrt((1 - exp(-2 * kappa_theta * dt)) / (2 * kappa_theta)) * (rho * dW_xi + sqrt(1 - rho^2) * dW_theta)
    #xim <- xim + kappa_xi * (thetam - xim) * dt + sqrt(dt) * sigma_xi * dW_xi
    xim <- xim * exp(-kappa_xi * dt) +
      thetam * (1 - exp(-kappa_xi * dt)) +
      sigma_xi * sqrt((1 - exp(-2 * kappa_xi * dt)) / (2 * kappa_xi)) * dW_xi
    
    ### Generate deterministic jump size
    J_D_Normal <- rnorm(1, mean = gamma_Q + xim, sd = omega) #rnorm(1, mean = gamma_Q + t(Gamma_Q) %*% X_j[, i - 1], sd = omega)
    
    if (any(t == jumps)) { #(any(abs(t - jumps) < 1e-5))
      ### Increment the counting process
      N_D <- 1
    } else {
      N_D <- 0
    }
    
    ### IORB Specific Short-Rate
    r[i] <- r[i - 1] + J_P_Normal + J_D_Normal * N_D
  }
  return(sum(r + s) * dt)
}

### Stochastic SOFR Spread + Normal Scheduled Jumps + Mixed-Exponential Unscheduled Jumps
SimNME <- function() {
  s <- numeric(n_steps + 1)
  s[1] <- s0
  r <- numeric(n_steps + 1)
  r[1] <- r0
  for (i in 2:(n_steps + 1)) {
    t <- i * dt
    #print(i)
    #s[i] <- s[i - 1] + kappa_s_jump * (theta_s_jump - s[i - 1]) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
    s[i] <- s[i - 1] * exp(-kappa_s_jump * dt) +
      theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
      sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
    
    ### Unscheduled jumps
    N_P <- rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * 1/360)
    ### Simulate jump sizes Z_1,Z_2,...
    if (N_P > 0) {
      U <- runif(N_P)                                              ### Uniform to determine component
      #Z_u <- rexp(N_P, rate = sum(lambda_J * eta_J))               ### Positive jump sizes (Exp(eta))
      #Z_d <- -rexp(N_P, rate = sum(q_J * (-theta_J)))              ### Negative jump sizes (Exp(theta)), negated
      #Z_u <- rexp(N_P, rate = eta_J)                               ### Positive jump sizes (Exp(eta))
      #Z_d <- -rexp(N_P, rate = -theta_J)                           ### Negative jump sizes (Exp(theta)), negated
      Jumpi <- sample(1:length(lambda_J), 1, prob = lambda_J)      ### Select component
      Z_u <- rexp(N_P, rate = eta_J[Jumpi])                        ### Positive jump sizes (Exp(eta))
      Jumpj <- sample(1:length(lambda_J), 1, prob = lambda_J)      ### Select component
      Z_d <- -rexp(N_P, rate = -theta_J[Jumpj])                    ### Negative jump sizes (Exp(theta)), negated
      
      ### Assign based on probability
      Z_MixedExp <- ifelse(U < pu_J, Z_u, Z_d)
    } else {
      Z_MixedExp <- 0
    }
    ### Update the jump process
    J_P_MixedExp <- sum(Z_MixedExp) ### Compute compound Poisson sum
    
    ### Scheduled jumps
    dW_xi <- rnorm(1)
    dW_theta <- rnorm(1)
    #thetam <- thetam + kappa_theta * (theta_theta - thetam) * dt + sqrt(dt) * sigma_theta * (rho * dW_xi + sqrt(1 - rho^2) * dW_theta)
    thetam <- thetam * exp(-kappa_theta * dt) +
      theta_theta * (1 - exp(-kappa_theta * dt)) +
      sigma_theta * sqrt((1 - exp(-2 * kappa_theta * dt)) / (2 * kappa_theta)) * (rho * dW_xi + sqrt(1 - rho^2) * dW_theta)
    #xim <- xim + kappa_xi * (thetam - xim) * dt + sqrt(dt) * sigma_xi * dW_xi
    xim <- xim * exp(-kappa_xi * dt) +
      thetam * (1 - exp(-kappa_xi * dt)) +
      sigma_xi * sqrt((1 - exp(-2 * kappa_xi * dt)) / (2 * kappa_xi)) * dW_xi
    
    ### Generate deterministic jump size
    J_D_Normal <- rnorm(1, mean = gamma_Q + xim, sd = omega) #rnorm(1, mean = gamma_Q + t(Gamma_Q) %*% X_j[, i - 1], sd = omega)
    
    if (any(t == jumps)) { #(any(abs(t - jumps) < 1e-5))
      ### Increment the counting process
      N_D <- 1
    } else {
      N_D <- 0
    }
    
    ### IORB Specific Short-Rate
    r[i] <- r[i - 1] + J_P_MixedExp + J_D_Normal * N_D
  }
  return(sum(r + s) * dt)
}

### Stochastic SOFR Spread + Skellam Scheduled Jumps + Normal Unscheduled Jumps
SimSN <- function() {
  s <- numeric(n_steps + 1)
  s[1] <- s0
  r <- numeric(n_steps + 1)
  r[1] <- r0
  for (i in 2:(n_steps + 1)) {
    t <- i * dt
    #print(i)
    #s[i] <- s[i - 1] + kappa_s_jump * (theta_s_jump - s[i - 1]) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
    s[i] <- s[i - 1] * exp(-kappa_s_jump * dt) +
      theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
      sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
    
    ### Unscheduled jumps
    N_P <- rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * 1/360)
    ### Simulate jump sizes Z_1,Z_2,...
    if (N_P > 0) {
      Z_Normal <- rnorm(N_P, mean = mu_P, sd = sigma_P)
    } else {
      Z_Normal <- 0
    }
    ### Update the jump process
    J_P_Normal <- sum(Z_Normal) ### Compute compound Poisson sum
    
    ### Scheduled jumps
    dW_zeta <- rnorm(1)
    zetam <- max(zetam + kappa_zeta * (theta_zeta - zetam) * dt + sigma_zeta * sqrt(pmax(zetam, 0)) * sqrt(dt) * dW_zeta +
                   1 / 4 * sigma_zeta^2 * dt * (dW_zeta^2 - 1), 0)
    
    ### Generate deterministic jump size
    J_D_Skellam <- nu_D * (rpois(1, lambda = gamma_Q + mu_d * zetam) - rpois(1, lambda = gamma_Q + mu_d * zetam))
    
    if (any(t == jumps)) { #(any(abs(t - jumps) < 1e-5))
      ### Increment the counting process
      N_D <- 1
    } else {
      N_D <- 0
    }
    
    ### IORB Specific Short-Rate
    r[i] <- r[i - 1] + J_P_Normal + J_D_Skellam * N_D
  }
  return(sum(r + s) * dt)
}

### Stochastic SOFR Spread + Skellam Scheduled Jumps + Mixed-Exponential Unscheduled Jumps
SimSME <- function() {
  s <- numeric(n_steps + 1)
  s[1] <- s0
  r <- numeric(n_steps + 1)
  r[1] <- r0
  for (i in 2:(n_steps + 1)) {
    t <- i * dt
    #print(i)
    #s[i] <- s[i - 1] + kappa_s_jump * (theta_s_jump - s[i - 1]) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
    s[i] <- s[i - 1] * exp(-kappa_s_jump * dt) +
      theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
      sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
    
    ### Unscheduled jumps
    N_P <- rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * dt) #rpois(1, lambda = nu_P * 1/360)
    ### Simulate jump sizes Z_1,Z_2,...
    if (N_P > 0) {
      U <- runif(N_P)                                              ### Uniform to determine component
      #Z_u <- rexp(N_P, rate = sum(lambda_J * eta_J))               ### Positive jump sizes (Exp(eta))
      #Z_d <- -rexp(N_P, rate = sum(q_J * (-theta_J)))              ### Negative jump sizes (Exp(theta)), negated
      #Z_u <- rexp(N_P, rate = eta_J)                               ### Positive jump sizes (Exp(eta))
      #Z_d <- -rexp(N_P, rate = -theta_J)                           ### Negative jump sizes (Exp(theta)), negated
      Jumpi <- sample(1:length(lambda_J), 1, prob = lambda_J)      ### Select component
      Z_u <- rexp(N_P, rate = eta_J[Jumpi])                        ### Positive jump sizes (Exp(eta))
      Jumpj <- sample(1:length(lambda_J), 1, prob = lambda_J)      ### Select component
      Z_d <- -rexp(N_P, rate = -theta_J[Jumpj])                    ### Negative jump sizes (Exp(theta)), negated
      
      ### Assign based on probability
      Z_MixedExp <- ifelse(U < pu_J, Z_u, Z_d)
    } else {
      Z_MixedExp <- 0
    }
    ### Update the jump process
    J_P_MixedExp <- sum(Z_MixedExp) ### Compute compound Poisson sum
    
    ### Scheduled jumps
    dW_zeta <- rnorm(1)
    zetam <- max(zetam + kappa_zeta * (theta_zeta - zetam) * dt + sigma_zeta * sqrt(pmax(zetam, 0)) * sqrt(dt) * dW_zeta +
                   1 / 4 * sigma_zeta^2 * dt * (dW_zeta^2 - 1), 0)
    
    ### Generate deterministic jump size
    J_D_Skellam <- nu_D * (rpois(1, lambda = gamma_Q + mu_d * zetam) - rpois(1, lambda = gamma_Q + mu_d * zetam))
    
    if (any(t == jumps)) { #(any(abs(t - jumps) < 1e-5))
      ### Increment the counting process
      N_D <- 1
    } else {
      N_D <- 0
    }
    
    ### IORB Specific Short-Rate
    r[i] <- r[i - 1] + J_P_MixedExp + J_D_Skellam * N_D
  }
  return(sum(r + s) * dt)
}




####################################################################################################################
####################################################################################################################
#---------------------------------- Calculations -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Run simulations
SimSOFR <- replicate(n_sim, SimSOFRSpread())
SimSOFRNN <- replicate(n_sim, SimNN())
SimSOFRSN <- replicate(n_sim, SimSN())
SimSOFRNME <- replicate(n_sim, SimNME())
SimSOFRSME <- replicate(n_sim, SimSME())



####################################################################################################################
####################################################################################################################
#------------------------------------- Plot ------------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Create data frame for plotting
SimPDFDF <- data.frame(
  IntegratedRate = c(SimSOFR, SimSOFRNN, SimSOFRNME, SimSOFRSN, SimSOFRSME),
  Model = factor(rep(c("Diffusion", "NN", "NME", "SN", "SME"),
                     each = n_sim))
)

### Plot
ggplot(SimPDFDF, aes(x = IntegratedRate, color = Model)) +
  geom_density(size = 1) +
  scale_color_manual(name = "Jump Distribution", values = c("NN" = "#901a1E", "NME" = "steelblue", "SN" = "#39641c", "SME" = "#666666", "Diffusion" = "hotpink3")) +
  labs(title = "Probability Density Function of Jump Model",
       x = "Integrated Interest Rate", y = "PDF") +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")

