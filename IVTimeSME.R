####################################################################################################################
####################################################################################################################
#------------------------------------ Settings ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Parameters
T <- 1/12
S <- 0
nsteps <- 10000
dt <- T / nsteps
Daily <- seq(1, T * 360) / 360  ### Daily steps

### Initial Values
zeta0 <- theta_zeta
s0 <- 0.003
r0 <- 0.03

### Jump
jump <- 0/360#10/360
jump_ <- 0/360#9/360


####################################################################################################################
####################################################################################################################
#--------------------------------- 1M SOFR Futures Price -----------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Expectation and Futures Price
SOFRSpot_MixedExpSkellam <- function(S, t, T, jump, jump_, zetat, rt, st, intS_t){
  if (t < jump){
    1/(T - S) * ((T - t) * rt +                                                       ### Constant Rate Contribution
                   1 / kappa_s_jump * (1 - exp(-kappa_s_jump * (T - t))) * st +       ### SOFR Spread Contribution
                   theta_s_jump /kappa_s_jump * (exp(-kappa_s_jump * (T - t)) - 1) +  ### SOFR Spread Contribution
                   theta_s_jump * (T - t) +                                           ### SOFR Spread Contribution
                   1/2 * (T - t)^2 * nu_P * (pu_J * sum(lambda_J * 1 / eta_J) -       ### Unscheduled Jumps Contribution
                                               qd_J * sum(q_J * 1 / theta_J)) +       ### Unscheduled Jumps Contribution
                   (T - jump) * nu_D *                                                                   ### Scheduled Jumps Contribution
                   #(xit_u * exp(-kappa_u * (jump_ - t)) + theta_u * (1 - exp(-kappa_u * (jump_ - t))) -  ### Scheduled Jumps Contribution
                   #   xit_d * exp(-kappa_d * (jump_ - t)) + theta_d * (1 - exp(-kappa_d * (jump_ - t)))) + ### Scheduled Jumps Contribution
                   (GammaQu * MeanXSkellam(t, jump_, zetat) -             ### Scheduled Jumps Contribution
                      GammaQd * MeanXSkellam(t, jump_, zetat)) +          ### Scheduled Jumps Contribution
                   intS_t)                                                            ### Integral from S to t
  } else if (t >= jump) {
    1/(T - S) * ((T - t) * rt +                                                       ### Constant Rate Contribution
                   1 / kappa_s_jump * (1 - exp(-kappa_s_jump * (T - t))) * st +       ### SOFR Spread Contribution
                   theta_s_jump /kappa_s_jump * (exp(-kappa_s_jump * (T - t)) - 1) +  ### SOFR Spread Contribution
                   theta_s_jump * (T - t) +                                           ### SOFR Spread Contribution
                   1/2 * (T - t)^2 * nu_P * (pu_J * sum(lambda_J * 1 / eta_J) -       ### Unscheduled Jumps Contribution
                                               qd_J * sum(q_J * 1 / theta_J)) +       ### Unscheduled Jumps Contribution
                   intS_t)                                                            ### Integral from S to t
    #1/(T - S) * (#(T - S) * r[(t*360) + 1] +
    #  theta_s_jump * (T - S) + 
    #    (exp(-kappa_s_jump * (S - t)) - exp(-kappa_s_jump * (T - t)))/kappa_s_jump * 
    #    (s - theta_s_jump))
  }      
}


### Strike
K_MixedExpSkellam <- 100 * (1 - SOFRSpot_MixedExpSkellam(S = 0, t = 0, T = T, jump = jump, jump_ = jump_, zetat = zeta0, rt = r0, st = s0, intS_t = 0))
k_MixedExpSkellam <- (T - S) * (100 - K_MixedExpSkellam) / 100



####################################################################################################################
####################################################################################################################
#------------------------------------- q_hat Function --------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Exact solution for beta_q
beta_q <- function(tau, x){
  (1 + x)/kappa_s_jump * (exp(-kappa_s_jump * tau) - 1)
}

### Exact solution for alpha_q
alpha_q <- function(tau, x){
  (1 + x) * theta_s_jump/kappa_s_jump * (1 - exp(-kappa_s_jump * tau)) -
    theta_s_jump * (1 + x) * tau +
    (sigma_s_jump^2 * (1 + x)^2)/(2 * kappa_s_jump^2) * tau +
    (sigma_s_jump^2 * (1 + x)^2)/(4 * kappa_s_jump^3) * (1 - exp(-2 * kappa_s_jump * tau)) -
    (sigma_s_jump^2 * (1 + x)^2)/(kappa_s_jump^3) * (1 - exp(-kappa_s_jump * tau))
}


alpha_qJ_MixedExp <- function(tau, x){
  #nu_P * (- tau +
  #          pu_J * sum((lambda_J * eta_J) / ((1 + x) + kappa_s_jump * eta_J) * log((((1 + x) + kappa_s_jump * eta_J) * exp(kappa_s_jump * tau) - (1 + x)) / (kappa_s_jump * eta_J))) -
  #          qd_J * sum((q_J * theta_J) / ((1 + x) - kappa_s_jump * theta_J) * log((((1 + x) - kappa_s_jump * theta_J) * exp(kappa_s_jump * tau) - (1 + x)) / (-kappa_s_jump * theta_J))))
  Pos <- numeric(length(lambda_J))
  for(i in 1:length(lambda_J)){
    #Pos[i] <- pu_J * (lambda_J[i] * eta_J[i]) / ((1 + x) + kappa_s_jump * eta_J[i]) *
    #  log((((1 + x) + kappa_s_jump * eta_J[i]) * exp(kappa_s_jump * tau) - (1 + x)) /
    #        (kappa_s_jump * eta_J[i]))
    Pos[i] <- pu_J * (lambda_J[i] * eta_J[i]) / (1 + x) *
      log(((1 + x) * tau + eta_J[i]) / (eta_J[i]))
  }
  Neg <- numeric(length(q_J))
  for(i in 1:length(q_J)){
    #Neg[i] <- qd_J * (q_J[i] * theta_J[i]) / ((1 + x) - kappa_s_jump * theta_J[i]) *
    #  log((((1 + x) - kappa_s_jump * theta_J[i]) * exp(kappa_s_jump * tau) - (1 + x)) /
    #        (-kappa_s_jump * theta_J[i]))
    Neg[i] <- qd_J * (q_J[i] * theta_J[i]) / (1 + x) *
      log((-(1 + x) * tau + theta_J[i]) / (theta_J[i]))
  }
  return(nu_P * (-tau + sum(Pos) - sum(Neg)))
}

#beta_qu <- function(x, tau){
#  u <- mu_u * (exp(nu_D * (-(T - jump) * (1 + x))) - 1)
#  B <- sigma_u^2 / (2 * kappa_u)
#  A <- 1/u - B
#  
#  return(1 / (A * exp(kappa_u * tau) + B))
#  #return(1 / (sigma_u^2 / (2 * kappa_u) * (1 - exp(kappa_u * tau)) + 1 / u * exp(kappa_u * tau)))
#}
#
#alpha_qu <- function(x, tau){
#  u <- mu_u * (exp(nu_D * (-(T - jump) * (1 + x))) - 1)
#  B <- sigma_u^2 / (2 * kappa_u)
#  A <- 1/u - B
#  
#  return(theta_u * 1 / B * (kappa_u * tau -
#                              log((A * exp(kappa_u * tau) + B) / (A + B))))
#}
#
#beta_qd <- function(x, tau){
#  u <- mu_d * (exp(-nu_D * (-(T - jump) * (1 + x))) - 1)
#  B <- sigma_d^2 / (2 * kappa_d)
#  A <- 1/u - B
#  
#  return(1 / (A * exp(kappa_d * tau) + B))
#  #return(1 / (sigma_d^2 / (2 * kappa_d) * (1 - exp(kappa_d * tau)) + 1 / u * exp(kappa_d * tau)))
#}
#
#alpha_qd <- function(x, tau){
#  u <- mu_d * (exp(-nu_D * (-(T - jump) * (1 + x))) - 1)
#  B <- sigma_d^2 / (2 * kappa_d)
#  A <- 1/u - B
#  
#  return(theta_d * 1 / B * (kappa_d * tau -
#                              log((A * exp(kappa_d * tau) + B) / (A + B))))
#}

beta_qD <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-(T - jump) * (1 + x))) - 1) + mu_d * (exp(-nu_D * (-(T - jump) * (1 + x))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(1 / (A * exp(kappa_zeta * tau) + B))
  #return(1 / (sigma_zeta^2 / (2 * kappa_zeta) * (1 - exp(kappa_zeta * tau)) + 1 / u * exp(kappa_zeta * tau)))
}

alpha_qD <- function(x, tau){
  u <- mu_u * (exp(nu_D * (-(T - jump) * (1 + x))) - 1) + mu_d * (exp(-nu_D * (-(T - jump) * (1 + x))) - 1)
  B <- sigma_zeta^2 / (2 * kappa_zeta)
  A <- 1/u - B
  
  return(theta_zeta * 1 / B * (kappa_zeta * tau -
                              log((A * exp(kappa_zeta * tau) + B) / (A + B))))
}

### Scaling
h <- 1.3

q_hat_MixedExpSkellam <- function(w, S, t, T, jump, jump_, K, ZCB, zetat, rt, st, intS_t){ 
  x <- h + complex(real = 0, imaginary = 1) * w
  
  #------------------ Strike Exponential -----------------------------------------------------------------------
  k <- (T - S)/100 * (100 - K)
  
  #k_hat <- k - trapz(seq(S, t, by = 1/360), SOFRSampleWeekend$SOFR[((S*360) + 1):((t*360) + 1)])
  #k_hat <- k - sum((head(SOFRSampleWeekend$SOFR[((S*360) + 1):((t*360) + 1)], -1) + tail(SOFRSampleWeekend$SOFR[((S*360) + 1):((t*360) + 1)], -1)) / 2 * (1/360)) ### Riemann Approxomation by mid-point rule
  k_hat <- k - intS_t
  
  ### Compute Strike Term
  StrikeTerm <- exp(x * k_hat)
  
  #------------------ ZCB Price -----------------------------------------------------------------------
  ### ODE system for alpha_J
  #func_J <- function(alpha_J, tau){
  #  nu_P * (exp(-mu_P * tau + (sigma_P * tau)^2/2) - 1)
  #}
  
  ### Solve for alpha_J(tau) using runge.kutta
  #taus <- seq(0, T - t, length.out = 100)  # Integration grid up to tau
  #alpha_J_RK <- runge.kutta(f = func_J, initial = 0, x = taus)
  
  ### Extract the final value of alpha_q
  #alpha_J_val <- tail(alpha_J_RK, 1)             ### w/ Unscheduled Jumps
  #alpha_J_val <- 0                               ### w/ No Unscheduled Jumps
  
  ### Compute ZCB price
  #ZCB <- exp(-(T - t) * rt + alpha_J_val + alpha_s(T - t) + beta_s(T - t) * st)
  
  #------------------ First Term -----------------------------------------------------------------------
  ### Compute First Term
  FirstTerm <- exp(-(T - t) * (1 + x) * rt)
  #FirstTerm <- 1
  
  #------------------ First Expectation -----------------------------------------------------------------------
  ### Compute First Expectation
  if (t < jump){
    ### MFG for non-central chi-squared
    #FirstExpectation <- 1 / (1 - (exp(nu_D * (-(T - jump) * (1 + x))) - 1) * 2 / kappa_u * sigma_u^2 *
    #                           (1 - exp(-kappa_u * (jump_ - t))))^((2 * kappa_u * theta_u) / (sigma_u^2)) *
    #  exp(((exp(nu_D * (-(T - jump) * (1 + x))) - 1) * xit_u * exp(-kappa_u * (jump_ - t))) /
    #        (1 - (exp(nu_D * (-(T - jump) * (1 + x))) - 1) * 2 / kappa_u * sigma_u^2 * (1 - exp(-kappa_u * (jump_ - t))))) *
    #  1 / (1 - (exp(-nu_D * (-(T - jump) * (1 + x))) - 1) * 2 / kappa_d * sigma_d^2 *
    #         (1 - exp(-kappa_d * (jump_ - t))))^((2 * kappa_d * theta_d) / (sigma_d^2)) *
    #  exp(((exp(-nu_D * (-(T - jump) * (1 + x))) - 1) * xit_d * exp(-kappa_d * (jump_ - t))) /
    #        (1 - (exp(-nu_D * (-(T - jump) * (1 + x))) - 1) * 2 / kappa_d * sigma_d^2 * (1 - exp(-kappa_d * (jump_ - t)))))
    ### ODEs
    #FirstExpectation <- exp(alpha_qu(x, tau = jump_ - t) + beta_qu(x, tau = jump_ - t) * xit_u +
    #                          alpha_qd(x, tau = jump_ - t) + beta_qd(x, tau = jump_ - t) * xit_d)
    FirstExpectation <- exp(alpha_qD(x, tau = jump_ - t) + beta_qD(x, tau = jump_ - t) * zetat)
  } else if (t >= jump) {
    FirstExpectation <- 1                               ### w/ No Scheduled Jumps
  }
  
  #------------------ Second Expectation -----------------------------------------------------------------------
  ### Compute Second Expectation
  SecondExpectation <- exp(alpha_q(T - t, x) + alpha_qJ_MixedExp(T - t, x) + beta_q(T - t, x) * st)          ### Jump-part of alpha^q
  #SecondExpectation <- exp(alpha_q(T - t, x) + beta_q(T - t, x) * st)          ### Jump-part of alpha^q
  
  #------------------ Final Price -----------------------------------------------------------------------
  ### Align dimensions
  StrikeTerm <- as.vector(StrikeTerm)
  #ZCB <- as.vector(ZCB)
  FirstTerm <- as.vector(FirstTerm)
  FirstExpectation <- as.vector(FirstExpectation)
  SecondExpectation <- as.vector(SecondExpectation)
  
  return(StrikeTerm * 1/ZCB * FirstTerm * FirstExpectation * SecondExpectation)
}


####################################################################################################################
####################################################################################################################
#------------------------------------- Call Price Function ---------------------------------------------------------
####################################################################################################################
####################################################################################################################
price1MSOFROption_MixedExpSkellam <- function(S, t, T, jump, jump_, K, ZCB, zetat, rt, st, intS_t) { 
  umax <- 50000
  integrand_q_hat <- function(w) {
    Re(q_hat_MixedExpSkellam(w, S = S, t = t, T = T, jump = jump, jump_ = jump_, K = K, ZCB = ZCB, zetat = zetat, rt = rt, st = st, intS_t = intS_t) / 
         ((h + complex(real = 0, imaginary = 1) * w)^2))
  }
  
  int <- (1 / pi) * integrate(integrand_q_hat, 0, umax)$value
  #pi <- 100/(T - S) * ZCB_s(t = t, T = T, rt = rt, st = st) * int
  pi <- 100/(T - S) * ZCB * int
  
  return(pi)
}

####################################################################################################################
####################################################################################################################
#---------------------------------- Implied Volatility -------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Black-Scholes Pricing
BlackScholesFormula <- function(spot, strike, t, T, r, div, D, sigma){
  d1 = (1/(sigma * sqrt(T - t))) * (log(spot/strike) + (r + 0.5 * sigma^2) * (T - t))
  d2 = d1 - sigma * sqrt(T - t)
  pi_BS = D * (spot * exp(-div * T) * pnorm(d1) - exp(-r * (T - t)) * strike * pnorm(d2))
  return(pi_BS)
}

### Implied Volatility
IVTest <- function(spot, strike, t, T, D, c) {
  Diff <- function(sigmaIV) BlackScholesFormula(spot, strike, t, T, r = 0, div = 0, D, sigmaIV) - c
  
  ### Handling boundary conditions
  if (c >= D * spot) {
    return(Inf)
  } else if (c <= D * max(spot - strike, 0)) {
    return(0)
  } else {
    ### Initial guess
    sigmaIV <- 0.02
    
    ### Determine search bounds
    if (Diff(sigmaIV) > 0) {
      b <- sigmaIV
      sigmaIV <- sigmaIV / 2
      while (Diff(sigmaIV) > 0) {
        sigmaIV <- sigmaIV / 2
      }
      a <- 0.99 * sigmaIV
    } else {
      a <- sigmaIV
      sigmaIV <- sigmaIV * 2
      while (Diff(sigmaIV) < 0) {
        sigmaIV <- sigmaIV * 2
      }
      b <- 1.01 * sigmaIV
    }
    
    ### Use uniroot to find the implied volatility
    IV <- tryCatch(uniroot(Diff, lower = a, upper = b)$root, error = function(e) NA)
    return(IV)
  }
}




####################################################################################################################
####################################################################################################################
#---------------------------------- Calculations -------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Simulation of Master Path
set.seed(3333333) #set.seed(NULL), 23, 333, 3333333
int_sm <- 0
zetam <- zeta0
sm <- s0
rm <- r0
J_P <- 0
J_D <- 0
N_D <- 0
IV_BS_Bach_MixedExpSkellam <- IV_BS_MixedExpSkellam <- numeric(nsteps)
S0_MixedExpSkellam <- numeric(nsteps)

for (i in 1:(nsteps - 10)) {  ### Avoid issues in the last steps
  print(i)
  t <- i * dt
  
  ### Simulate sm at each fine time step
  #sm <- sm + kappa_s_jump * (theta_s_jump - sm) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
  sm <- sm * exp(-kappa_s_jump * dt) +
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
    Z_d <- -rexp(N_P, rate = theta_J[Jumpj])                    ### Negative jump sizes (Exp(theta)), negated
    
    ### Assign based on probability
    Z <- ifelse(U < pu_J, Z_u, Z_d)
  } else {
    Z <- 0
  }
  ### Update the jump process
  J_P <- sum(Z) ### Compute compound Poisson sum
  print(J_P)
  
  ### Scheduled jumps
  dW_zeta <- rnorm(1)
  #dW_theta_u <- rnorm(1)
  zetam <- max(zetam + kappa_zeta * (theta_zeta - zetam) * dt + sigma_zeta * sqrt(pmax(zetam, 0)) * sqrt(dt) * dW_zeta +
                 1 / 4 * sigma_zeta^2 * dt * (dW_zeta^2 - 1), 0)
  dW_xi_d <- rnorm(1)
  #dW_theta_d <- rnorm(1)
  #xim_d <- max(xim_d + kappa_d * (theta_d - xim_d) * dt + sigma_d * sqrt(pmax(xim_d, 0)) * sqrt(dt) * dW_xi_d +
  #               1 / 4 * sigma_d^2 * dt * (dW_xi_d^2 - 1), 0)
  #xim_d <- xim_u
  
  ### Generate deterministic jump size
  #J_D <- nu_D * (rpois(1, lambda = mu_u) - rpois(1, lambda = mu_d))
  J_D <- nu_D * (rpois(1, lambda = gamma_Q + mu_u * zetam) - rpois(1, lambda = gamma_Q + mu_d * zetam))
  
  if (t == jump) { #(time_series[i] %in% jump_dates) #(is_close(t, jump_dates))
    ### Increment the counting process
    N_D <- 1
  } else {
    N_D <- 0
  }
  print(J_D * N_D)
  
  ### Update rm and int_sm only at full-day steps
  if (any(abs(t - Daily) < 1e-6)) {  ### Check if t matches a full day
    ### Determine rt from SOFR r
    idx <- max(which(Daily <= t), na.rm = TRUE)  ### Find the last time point before or equal to t
    if (is.infinite(idx)) idx <- 1  ### Avoid empty selection
    #rm <- rm + J_P + J_D * N_D ### Update once a day
    
    ### Integral of r^s from S to t daily
    #int_sm <- int_sm + (rm + sm) * 1/360  ### Accumulate integral only at daily steps
  }
  rm <- rm + J_P + J_D * N_D ### Updating at every time step
  
  ### Integral of r^s from S to t for each t
  int_sm <- int_sm + (rm + sm) * dt ### Accumulate integral at every time step
  
  ### Call Price
  Call <- price1MSOFROption_MixedExpSkellam(S, t, T, jump, jump_, K_MixedExpSkellam, ZCB_s_SME(t, T, zetat = zetam, rt = rm, sm), zetat = zetam, rt = rm, st = sm, intS_t = int_sm)
  
  ### Spots
  S0_MixedExpSkellam[i] <- 100 * (1 - SOFRSpot_MixedExpSkellam(S, t, T, jump, jump_, zetat = zetam, rt = rm, st = sm, intS_t = int_sm))
  
  ### Implicit Volatility
  IV_BS_MixedExpSkellam[i] <- IVTest(S0_MixedExpSkellam[i], K_MixedExpSkellam, t, T, ZCB_s_SME(t, T, zetam, rm, sm), Call)
}

plot(seq(0, 1/12, length.out = nsteps), IV_BS_MixedExp, type = "l", col = "#901a1E", lwd = 1, ylim = c(0, 0.08),
     xlab = "t", ylab = "Implied Volatility", main = "Black-76 Implied Volatility for t in [0, 1/12]")
lines(seq(0, 1/12, length.out = nsteps), IV_BS_MixedExpSkellam, col = "#666666", lwd = 1)
abline(v = (10 - 0)/360, lty = 2, lwd = 1)
legend("topright",
       legend = c("Normal","Skellam", "Jump"),
       col = c("#901a1E", "#666666", "black"),
       title = expression("Dist. of " * J^D),
       lwd = c(1, 1, 1), lty = c(1, 1, 2), cex = 0.7, inset = 0.00)



### Plotting the results
plot(S0_MixedExpSkellam[1:(nsteps - 100)], col = "steelblue")
plot(seq(0, T, length.out = nsteps), IV_BS_Bach_MixedExpSkellam, type = "l", col = "steelblue", 
     xlab = "t", ylab = "Implied Volatility", main = "Imp Black Vol for t in [0, 1/12]")
abline(v = (10 - 0)/360, lty = 2)
plot(seq(0, T, length.out = nsteps), IV_BS_MixedExpSkellam, type = "l", col = "steelblue", 
     xlab = "t", ylab = "Implied Volatility", main = "Imp Black Vol for t in [0, 1/12]")
abline(v = (10 - 0)/360, lty = 2)
plot(seq(0, T, length.out = nsteps), IV_BS_MixedExpSkellam, type = "l", col = "steelblue", 
     xlab = "t", ylab = "Implied Volatility", main = "Imp Black Vol for t in [0, 1/12]", ylim = c(0, 0.01))
abline(v = (10 - 0)/360, lty = 2)





T <- 1/12
nsteps <- 10000

plot(seq(0, 1/12, length.out = nsteps), IV_BS_MixedExp, type = "l", col = "#901a1E", lwd = 1, ylim = c(0, 0.07),
     xlab = "t", ylab = "Implied Volatility", main = "Black-76 Implied Volatility for t in [0, 1/12]")
lines(seq(0, 1/12, length.out = nsteps), IV_BS_Normal, col = "steelblue", lwd = 1)
lines(seq(0, 1/12, length.out = nsteps), IV_BS_MixedExpSkellam, col = "#666666", lwd = 1)
lines(seq(0, 1/12, length.out = nsteps), IV_BS_NormalSkellam, col = "#39641c", lwd = 1)
abline(v = (10 - 0)/360, lty = 2, lwd = 1)
legend("topright",
       legend = expression(
         paste(italic(N), "/", italic(ME)),
         paste(italic(N), "/", italic(N)),
         paste(italic(S), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         "Jump"),
       col = c("#901a1E", "steelblue", "#666666", "#39641c", "black"),
       title = expression("Dist. of " * J^D * " / " * J^P),
       lwd = c(1, 1, 1, 1, 1), lty = c(1, 1, 1, 1, 2), cex = 0.7, inset = 0.00)
legend("topright",
       legend = c("Normal/Mixed-Exp", "Normal/Normal", "Skellam/Mixed-Exp", "Skellam/Normal", "Scheduled Jump"),
       col = c("#901a1E", "steelblue", "#666666", "#39641c", "black"), title = "Dist. of J_D / J_P",
       lwd = c(1, 1, 1, 1, 1), lty = c(1, 1, 1, 1, 2), cex=0.7, inset = 0.00)



