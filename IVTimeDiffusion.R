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
s0 <- 0.003
r0 <- 0.03

### Jump
jump <- 10/360
jump_ <- 9/360


####################################################################################################################
####################################################################################################################
#--------------------------------- 1M SOFR Futures Price -----------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Expectation and Futures Price
SOFRSpot_Diffusion <- function(S, t, T, jump, jump_, rt, st, intS_t){
  if (t < jump){
    1/(T - S) * ((T - t) * rt +                                                       ### Constant Rate Contribution
                   1 / kappa_s_jump * (1 - exp(-kappa_s_jump * (T - t))) * st +       ### SOFR Spread Contribution
                   theta_s_jump /kappa_s_jump * (exp(-kappa_s_jump * (T - t)) - 1) +  ### SOFR Spread Contribution
                   theta_s_jump * (T - t) +                                           ### SOFR Spread Contribution
                   intS_t)                                                            ### Integral from S to t
  } else if (t >= jump) {
    1/(T - S) * ((T - t) * rt +                                                       ### Constant Rate Contribution
                   1 / kappa_s_jump * (1 - exp(-kappa_s_jump * (T - t))) * st +       ### SOFR Spread Contribution
                   theta_s_jump /kappa_s_jump * (exp(-kappa_s_jump * (T - t)) - 1) +  ### SOFR Spread Contribution
                   theta_s_jump * (T - t) +                                           ### SOFR Spread Contribution
                   intS_t)                                                            ### Integral from S to t
    #1/(T - S) * (#(T - S) * r[(t*360) + 1] +
    #  theta_s_jump * (T - S) + 
    #    (exp(-kappa_s_jump * (S - t)) - exp(-kappa_s_jump * (T - t)))/kappa_s_jump * 
    #    (s - theta_s_jump))
  }      
}

### Strike
K_Diffusion <- 100 * (1 - SOFRSpot_Diffusion(S = 0, t = 0, T = T, jump = jump, jump_ = jump_, rt = r0, st = s0, intS_t = 0))
k_Diffusion <- (T - S) * (100 - K_Diffusion) / 100


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

### Scaling
h <- 1.3

q_hat_Diffusion <- function(w, S, t, T, jump, jump_, K, ZCB, rt, st, intS_t){ 
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
    FirstExpectation <- 1
  } else if (t >= jump) {
    FirstExpectation <- 1                               ### w/ No Scheduled Jumps
  }
  
  #------------------ Second Expectation -----------------------------------------------------------------------
  ### Compute Second Expectation
  SecondExpectation <- exp(alpha_q(T - t, x) + beta_q(T - t, x) * st)          ### Jump-part of alpha^q
  
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
price1MSOFROption_Diffusion <- function(S, t, T, jump, jump_, K, ZCB, rt, st, intS_t) { 
  umax <- 50000
  integrand_q_hat <- function(w) {
    Re(q_hat_Diffusion(w, S = S, t = t, T = T, jump = jump, jump_ = jump_, K = K, ZCB = ZCB, rt = rt, st = st, intS_t = intS_t) / 
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
    #IV <- tryCatch(uniroot(Diff, lower = 1e-1113, upper = 1e-13)$root, error = function(e) NA)
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
sm <- s0
rm <- r0
J_P <- 0
J_D <- 0
N_D <- 0
IV_BS_Diffusion <- numeric(nsteps)
S0_Diffusion <- numeric(nsteps)

for (i in 1:(nsteps - 10)) {  ### Avoid issues in the last steps
  print(i)
  t <- i * dt
  
  ### Simulate sm at each fine time step
  #sm <- sm + kappa_s_jump * (theta_s_jump - sm) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
  sm <- sm * exp(-kappa_s_jump * dt) +
    theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
    sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(1)
  
  ### Unscheduled jumps
  J_P <- 0
  print(J_P)
  
  ### Scheduled jumps
  J_D <- 0
  
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
  Call <- price1MSOFROption_Diffusion(S, t, T, jump, jump_, K_Diffusion, ZCB_s(t, T, rt = rm, sm), rt = rm, st = sm, intS_t = int_sm)
  
  ### Spots
  S0_Diffusion[i] <- 100 * (1 - SOFRSpot_Diffusion(S, t, T, jump, jump_, rt = rm, st = sm, intS_t = int_sm))
  
  ### Implicit Volatility
  IV_BS_Diffusion[i] <- IVTest(S0_Diffusion[i], K_Diffusion, t, T, ZCB_s(t, T, rm, sm), Call)
}





