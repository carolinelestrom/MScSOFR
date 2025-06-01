####################################################################################################################
####################################################################################################################
#------------------------------------ Settings ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Parameters
T <- 12/12
S <- 11/12
nsteps <- 10000
dt <- T / nsteps

### Initial Values
zeta0 <- theta_zeta
s0 <- 0.003
r0 <- 0.03

### Jumps for 1Y (2025)
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
#--------------------------------- 1M SOFR Futures Price -----------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Expectation and Futures Price
SOFRSpot_Smile_NormalSkellam <- function(S, t, T, jumps, jumps_, zetat, rt, st, intS_t){
  1/(T - S) * (rt * (T - S) +                                                         ### Constant Rate Contribution
                 theta_s_jump * (T - S) +                                             ### SOFR Spread Contribution
                 (exp(-kappa_s_jump * (S - t)) - exp(-kappa_s_jump * (T - t))) /      ### SOFR Spread Contribution
                 kappa_s_jump * (st - theta_s_jump) +                                 ### SOFR Spread Contribution
                 1/2 * (T - S)^2 * nu_P * (pu_J * sum(lambda_J * 1 / eta_J) -         ### Unscheduled Jumps Contribution
                                             qd_J * sum(q_J * 1 / theta_J)) +         ### Unscheduled Jumps Contribution
                 min(T - S, T - jumps[1]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[1], zetat) -          ### 1st Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[1], zetat)) +       ### 1st Scheduled Jump
                 min(T - S, T - jumps[2]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[2], zetat) -          ### 2nd Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[2], zetat)) +       ### 2nd Scheduled Jump
                 min(T - S, T - jumps[3]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[3], zetat) -          ### 3rd Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[3], zetat)) +       ### 3rd Scheduled Jump
                 min(T - S, T - jumps[4]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[4], zetat) -          ### 4th Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[4], zetat)) +       ### 4th Scheduled Jump
                 min(T - S, T - jumps[5]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[5], zetat) -          ### 5th Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[5], zetat)) +       ### 5th Scheduled Jump
                 min(T - S, T - jumps[6]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[6], zetat) -          ### 6th Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[6], zetat)) +       ### 6th Scheduled Jump
                 min(T - S, T - jumps[7]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[7], zetat) -          ### 7th Scheduled Jump
                                                      GammaQd * MeanXSkellam(t, jumps_[7], zetat)) +       ### 7th Scheduled Jump
                 min(T - S, T - jumps[8]) * nu_D * (GammaQu * MeanXSkellam(t, jumps_[8], zetat) -          ### 8th Scheduled Jump in reference period
                                                      GammaQd * MeanXSkellam(t, jumps_[8], zetat)) +       ### 8th Scheduled Jump in reference period
                 intS_t)                                                              ### Integral from S to t
}


### Strike
K_Smile_NormalSkellam <- 100 * (1 - SOFRSpot_Smile_NormalSkellam(S = S, t = 0, T = T, jumps = jumps, jumps_ = jumps_, zetat = zeta0, rt = r0, st = s0, intS_t = 0))
k_Smile_NormalSkellam <- (T - S) * (100 - K_Smile_NormalSkellam) / 100

####################################################################################################################
####################################################################################################################
#------------------------------------- q_hat Function --------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Exact solution for beta_qs
beta_qs_ES <- function(x, tau){
  (1 + x) / kappa_s_jump * (exp(-kappa_s_jump * tau) - 1)
}

### Exact solution for \tilde{beta_qs}
beta_tilde_qs_EtES <- function(x, tau, S, T){
  1 / kappa_s_jump * (exp(-kappa_s_jump * tau) - 1) + beta_qs_ES(x, T - S) * exp(-kappa_s_jump * tau)
}

### Exact solution for alpha_qs
alpha_qs_ES <- function(x, tau){
  - theta_s_jump * (1 + x) * tau + ((sigma_s_jump^2 * (1 + x)^2) / (2 * kappa_s_jump^2)) * tau +
    ((1 + x) * theta_s_jump / kappa_s_jump - (sigma_s_jump^2 * (1 + x)^2) / (kappa_s_jump^3)) * (1 - exp(-kappa_s_jump * tau)) +
    ((sigma_s_jump^2 * (1 + x)^2) / (4 * kappa_s_jump^3)) * (1 - exp(-2 * kappa_s_jump * tau))
}

### Exact solution for \tilde{alpha_qs}
alpha_tilde_qs_EtES <- function(x, tau, S, T){
  #alpha_qs_ES(x, T - S) +
  - theta_s_jump * tau + ((sigma_s_jump^2) / (2 * kappa_s_jump^2)) * tau +
    (theta_s_jump / kappa_s_jump + theta_s_jump * beta_qs_ES(x, T - S) -
       (sigma_s_jump^2 / kappa_s_jump^3) - (sigma_s_jump^2 / kappa_s_jump^2) * beta_qs_ES(x, T - S)) * (1 - exp(-kappa_s_jump * tau)) +
    ((sigma_s_jump^2 / (4 * kappa_s_jump^3)) + (sigma_s_jump^2 / (2 * kappa_s_jump^2)) * beta_qs_ES(x, T - S) +
       (sigma_s_jump^2 / (4 * kappa_s_jump)) * (beta_qs_ES(x, T - S))^2) * (1 - exp(-2 * kappa_s_jump * tau))
}

### Scaling
h <- 1.3

q_hat_Smile_NormalSkellam <- function(w, S, t, T, jumps, jumps_, K, ZCB, zetat, rt, st, intS_t){ 
  x <- h + complex(real = 0, imaginary = 1) * w
  
  #------------------ Strike Exponential -----------------------------------------------------------------------
  k <- (T - S)/100 * (100 - K)
  
  ### Compute Strike Term
  StrikeTerm <- exp(x * k)
  
  #------------------ First Term -----------------------------------------------------------------------
  ### Compute First Term
  FirstTerm <- exp(-(S - t) * rt) * exp(-(T - S) * (1 + x) * rt)
  
  #------------------ First Expectation -----------------------------------------------------------------------
  ### Compute First Expectation
  FirstExpectation <- exp(alpha_q_tau8(x, tau = jumps_[8] - jumps_[7]) + 
                            alpha_q_tau7(x, tau = jumps_[7] - jumps_[6]) + 
                            alpha_q_tau6(x, tau = jumps_[6] - jumps_[5]) + 
                            alpha_q_tau5(x, tau = jumps_[5] - jumps_[4]) + 
                            alpha_q_tau4(x, tau = jumps_[4] - jumps_[3]) + 
                            alpha_q_tau3(x, tau = jumps_[3] - jumps_[2]) + 
                            alpha_q_tau2(x, tau = jumps_[2] - jumps_[1]) + 
                            alpha_q_tau1(x, tau = jumps_[1] - t) + 
                            beta_q_tau1(x, tau = jumps_[1] - t) * zetat)
  
  #------------------ Second Expectation -----------------------------------------------------------------------
  ### Define func_qJ_ES to compute alpha_qJ using the ODE
  func_qJ_ES <- function(tau, alpha_qJ_ES) {
    nu_P * (exp(mu_P * (- (1 + x) * tau) + (sigma_P * (- (1 + x) * tau))^2 / 2) - 1)
  }
  
  ### Use a simple RK4 or your own integration routine (e.g., from pracma or deSolve)
  taus_ES <- seq(0, T - S, length.out = 100)
  alpha_qJ_ES_vals <- numeric(length(taus_ES))
  alpha_qJ_ES_vals[1] <- 0                       ### Initial value
  dtau_ES <- taus_ES[2] - taus_ES[1]                  ### Time step difference
  
  for (i in 2:length(taus_ES)) {
    #k1_ES <- dtau_ES * func_qJ_ES(taus_ES[i - 1], alpha_qJ_ES_vals[i - 1])
    #k2_ES <- dtau_ES * func_qJ_ES(taus_ES[i - 1] + dtau_ES/2, alpha_qJ_ES_vals[i - 1] + k1_ES / 2)
    #k3_ES <- dtau_ES * func_qJ_ES(taus_ES[i - 1] + dtau_ES/2, alpha_qJ_ES_vals[i - 1] + k2_ES / 2)
    #k4_ES <- dtau_ES * func_qJ_ES(taus_ES[i - 1] + dtau_ES, alpha_qJ_ES_vals[i - 1] + k3_ES)
    #alpha_qJ_ES_vals[i] <- alpha_qJ_ES_vals[i - 1] + (k1_ES + 2*k2_ES + 2*k3_ES + k4_ES) / 6
    k1_ES <- func_qJ_ES(taus_ES[i - 1], alpha_qJ_ES_vals[i - 1])
    k2_ES <- func_qJ_ES(taus_ES[i - 1] + dtau_ES/2, alpha_qJ_ES_vals[i - 1] + dtau_ES * k1_ES / 2)
    k3_ES <- func_qJ_ES(taus_ES[i - 1] + dtau_ES/2, alpha_qJ_ES_vals[i - 1] + dtau_ES * k2_ES / 2)
    k4_ES <- func_qJ_ES(taus_ES[i - 1] + dtau_ES, alpha_qJ_ES_vals[i - 1] + dtau_ES * k3_ES)
    alpha_qJ_ES_vals[i] <- alpha_qJ_ES_vals[i - 1] + (dtau_ES / 6) * (k1_ES + 2*k2_ES + 2*k3_ES + k4_ES)
  }
  
  alpha_qJ_ES_val <- tail(alpha_qJ_ES_vals, 1)             ### Unscheduled Jump Contribution
  
  ### Define func_tilde_qJ_EtES to compute \tilde{alpha}_qJ using the ODE
  func_tilde_qJ_EtES <- function(tau, alpha_tilde_qJ_EtES) {
    nu_P * (exp(mu_P * ((1 + x) * (T - S) - tau) + (sigma_P * ((1 + x) * (T - S) - tau))^2 / 2) - 1)
  }
  
  ### Use a simple RK4 or your own integration routine (e.g., from pracma or deSolve)
  taus_EtES <- seq(0, S - t, length.out = 100)
  alpha_tilde_qJ_EtES_vals <- numeric(length(taus_EtES))
  alpha_tilde_qJ_EtES_vals[1] <- 0                       ### Initial value
  dtau_EtES <- taus_EtES[2] - taus_EtES[1]                  ### Time step difference
  
  for (i in 2:length(taus_EtES)) {
    #k1_EtES <- dtau_EtES * func_tilde_qJ_EtES(taus_EtES[i - 1], alpha_tilde_qJ_EtES_vals[i - 1])
    #k2_EtES <- dtau_EtES * func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES/2, alpha_tilde_qJ_EtES_vals[i - 1] + k1_EtES / 2)
    #k3_EtES <- dtau_EtES * func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES/2, alpha_tilde_qJ_EtES_vals[i - 1] + k2_EtES / 2)
    #k4_EtES <- dtau_EtES * func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES, alpha_tilde_qJ_EtES_vals[i - 1] + k3_EtES)
    #alpha_tilde_qJ_EtES_vals[i] <- alpha_tilde_qJ_EtES_vals[i - 1] + (k1_EtES + 2*k2_EtES + 2*k3_EtES + k4_EtES) / 6
    k1_EtES <- func_tilde_qJ_EtES(taus_EtES[i - 1], alpha_tilde_qJ_EtES_vals[i - 1])
    k2_EtES <- func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES/2, alpha_tilde_qJ_EtES_vals[i - 1] + dtau_EtES * k1_EtES / 2)
    k3_EtES <- func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES/2, alpha_tilde_qJ_EtES_vals[i - 1] + dtau_EtES * k2_EtES / 2)
    k4_EtES <- func_tilde_qJ_EtES(taus_EtES[i - 1] + dtau_EtES, alpha_tilde_qJ_EtES_vals[i - 1] + dtau_EtES * k3_EtES)
    alpha_tilde_qJ_EtES_vals[i] <- alpha_tilde_qJ_EtES_vals[i - 1] + (dtau_EtES / 6) * (k1_EtES + 2*k2_EtES + 2*k3_EtES + k4_EtES)
  }
  
  alpha_tilde_qJ_EtES_val <- tail(alpha_tilde_qJ_EtES_vals, 1)             ### Unscheduled Jump Contribution
  
  ### Compute Second Expectation
  SecondExpectation <- exp(alpha_qs_ES(x, T - S) + alpha_qJ_ES_val + ### initial for alpha
                             alpha_tilde_qs_EtES(x, tau = S - t, T = T, S = S) +
                             alpha_tilde_qJ_EtES_val +
                             beta_tilde_qs_EtES(x, tau = S - t, T = T, S = S) * st)
  
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
price1MSOFROption_Smile_NormalSkellam <- function(S, t, T, jumps, jumps_, K, ZCB, zetat, rt, st, intS_t) { 
  umax <- 50000
  integrand_q_hat <- function(w) {
    Re(q_hat_Smile_NormalSkellam(w, S = S, t = t, T = T, jumps = jumps, jumps_ = jumps_, K = K, ZCB = ZCB, zetat = zetat, rt = rt, st = st, intS_t = intS_t) / 
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
#---------------------------------- Calculations -------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
StrikeSmile <- seq(96-23, 96+23, 1) ### Strikes for the smile
CallSmile_NormalSkellam <- rep(NA, length(StrikeSmile))
IVSmile_NormalSkellam <- rep(NA, length(StrikeSmile))
for (i in 1:length(StrikeSmile)) {
  print(i)
  CallSmile_NormalSkellam[i] <- price1MSOFROption_Smile_NormalSkellam(S = S, t = 0, T = T, jumps = jumps, jumps_ = jumps_,
                                                                      K = StrikeSmile[i], ZCB_s_SN(t = 0, T = T, zetat = zeta0, rt = r0, st = s0),
                                                                      zetat = zeta0, rt = r0, st = s0, intS_t = 0)
  IVSmile_NormalSkellam[i] = IVTest(spot = K_Smile_NormalSkellam, strike = StrikeSmile[i],
                                    t = 0, T = T, ZCB_s_SN(t = 0, T = T, zetat = zeta0, rt = r0, st = s0), CallSmile_NormalSkellam[i])
}





### PLOT
### ====

plot(StrikeSmile, CallSmile_Normal, type = "l", lwd = 3, lty = 2,
     ylab = "Call Price", xlab = "Strike")
lines(StrikeSmile, CallSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
lines(StrikeSmile, CallSmile_MixedExpSkellam, col = "#666666", lwd = 3, lty = 2)
lines(StrikeSmile, CallSmile_NormalSkellam, col = "#39641c", lwd = 3, lty = 2)
lines(StrikeSmile, CallSmile, col = "hotpink3", lwd = 3, lty = 2)
lines(StrikeSmile, CallSmile_Test, col = "orchid4", lwd = 3, lty = 2)
lines(StrikeSmile, MCPrice, col = "#901a1E", lwd = 3, lty = 3)
abline(h = 0)
legend("topright",
       legend = expression(
         paste(italic(N), "/", italic(N)),
         paste(italic(N), "/", italic(ME)),
         paste(italic(S), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         "Diffusion",
         paste(italic(N), "/", italic(ME), " Fix Var"),
         "MC"
       ),
       col = c("black", "steelblue", "#666666", "#39641c", "hotpink3", "orchid4", "#901a1E"),
       title = "Call Price",
       lwd = c(1, 1, 1, 1, 1), lty = c(2, 2, 2, 2, 3), cex = 0.7, inset = 0.00)

plot(StrikeSmile, IVSmile_Normal, type = "l", lwd = 3, ylim = c(0, 0.3), lty = 2,
     ylab = "Black-76 Implicit Volatility", xlab = "Strike")
lines(StrikeSmile, IVSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
lines(StrikeSmile, IVSmile_MixedExpSkellam, col = "#666666", lwd = 3, lty = 2)
lines(StrikeSmile, IVSmile_NormalSkellam, col = "#39641c", lwd = 3, lty = 2)
lines(StrikeSmile, IVSmile, col = "hotpink3", lwd = 3, lty = 2)
lines(StrikeSmile, IVSmile_Test, col = "orchid4", lwd = 3, lty = 2)
lines(StrikeSmile, IVSmile_MC, col = "#901a1E", lwd = 3, lty = 3)
abline(h = 0)
legend("topright",
       legend = expression(
         paste(italic(N), "/", italic(N)),
         paste(italic(N), "/", italic(ME)),
         paste(italic(S), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         "Diffusion",
         paste(italic(N), "/", italic(ME), " Fix Var"),
         "MC"
       ),
       col = c("black", "steelblue", "#666666", "#39641c", "hotpink3", "orchid4", "#901a1E"),
       title = "Imp Vol",
       lwd = c(1, 1, 1, 1, 1, 1, 1), lty = c(2, 2, 2, 2, 2, 2, 3), cex = 0.7, inset = 0.00)








plot(StrikeSmile, CallSmile_Normal, type = "l", lwd = 3,
     ylab = "Call Price", xlab = "Strike")
lines(StrikeSmile, CallSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
abline(h = 0)
legend("topright", legend = c("Mixed-Exp", "Normal"),
       col = c("steelblue", "black", "#901a1E"), title = "Call Price",
       lwd = c(1, 1, 1), lty = c(2, 1, 3), cex=0.7, inset = 0.00)

plot(StrikeSmile, IVSmile_Normal, type = "l", lwd = 3,
     ylab = "Black-76 Implicit Volatility", xlab = "Strike")
lines(StrikeSmile, IVSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
abline(h = 0)
legend("bottomleft", legend = c("Mixed-Exp", "Normal"),
       col = c("steelblue", "black", "#901a1E"), title = "Imp Vol",
       lwd = c(1, 1, 1), lty = c(2, 1, 3), cex=0.7, inset = 0.00)




SmileData_Normal <- as.data.frame(cbind(StrikeSmile, CallSmile_Normal, IVSmile_Normal))
plot(StrikeSmile, IVSmile_Normal, type = "p")
plot(StrikeSmile, CallSmile_Normal, type = "l", lwd = 3)
lines(StrikeSmile, CallSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
lines(StrikeSmile, MCPrice, col = "#901a1E", lwd = 3, lty = 3)
abline(h = 0)
legend("topright", legend = c("Mixed-Exp", "Normal", "MC"),
       col = c("steelblue", "black", "#901a1E"), title = "Call Price",
       lwd = c(1, 1, 1), lty = c(2, 1, 3), cex=0.7, inset = 0.00)


plot(StrikeSmile, IVSmile, type = "l", col = "black")
lines(StrikeSmile, IVSmile_MixedExp, col = "#901a1E")
lines(StrikeSmile, IVSmile_Normal, col = "steelblue")

SmileData_Normal %>%
  ggplot(aes(x = StrikeSmile, y = IVSmile_Normal)) +
  geom_line(color = "#901a1E", size = 1.7) +
  theme_bw() +
  xlab("Strike") + ylab("Implied volatilities") +
  ggtitle("Implied volatility across strikes in the Bachelier model") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=23)) +
  theme(axis.text.x = element_text(size=17, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(50, 75, 100, 125, 150), limits = c(50, 150), 
                     labels = c(50, 75, "100 (ATM)", 125, 150)) +
  scale_y_continuous(breaks=c(0.12, 0.14, 0.16, 0.18, 0.20), limits = c(0.117, 0.21)) +
  geom_segment(aes(x = 100, y = 0.117, xend = 100, yend = 0.15), 
               linetype = "dotdash", linewidth = 1.3, color = "steelblue")  +
  geom_segment(aes(x = 50, y = 0.15, xend = 100, yend = 0.15), 
               linetype = "dotdash", linewidth = 1.3, color = "steelblue")



SmileData <- as.data.frame(cbind(StrikeSmile, CallSmile_MixedExp, IVSmile_MixedExp,
                                 CallSmile_Normal, IVSmile_Normal))
SmileData %>%
  ggplot(aes(x = StrikeSmile, y = IVSmile_MixedExp)) +
  geom_line(color = "#901a1E", size = 1.7) +
  geom_line(color = "steelblue", aes(y = IVSmile_Normal), size = 1.7) +
  theme_bw() +
  xlab("Strike") + ylab("Implied Volatility") +
  ggtitle("Black-76 Implied Volatility Across Strikes") +
  theme(plot.title = element_text(size=27, hjust=0)) +
  theme(axis.title = element_text(size=23)) +
  theme(axis.text.x = element_text(size=17, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=17)) +
  scale_x_continuous(breaks=c(96-23, 96-23/2, 96, 96+23/2, 96+23), limits = c(96-23, 96+23), 
                     labels = c(96-23, 96-23/2, "96 (~ATM)", 96+23/2, 96+23)) +
  scale_y_continuous(breaks=c(0.00, 0.05, 0.10, 0.15), limits = c(0.00, 0.15))


