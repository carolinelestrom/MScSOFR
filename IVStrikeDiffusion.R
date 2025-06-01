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
SOFRSpot_Smile <- function(S, t, T, rt, st, intS_t){
  1/(T - S) * (rt * (T - S) +                                                         ### Constant Rate Contribution
                 theta_s_jump * (T - S) +                                             ### SOFR Spread Contribution
                 (exp(-kappa_s_jump * (S - t)) - exp(-kappa_s_jump * (T - t))) /      ### SOFR Spread Contribution
                 kappa_s_jump * (st - theta_s_jump) +                                 ### SOFR Spread Contribution
                 intS_t)                                                              ### Integral from S to t
}


### Strike
K_Smile <- 100 * (1 - SOFRSpot_Smile(S = S, t = 0, T = T, rt = r0, st = s0, intS_t = 0))
k_Smile <- (T - S) * (100 - K_Smile) / 100

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

q_hat_Smile <- function(w, S, t, T, K, ZCB, rt, st, intS_t){ 
  x <- h + complex(real = 0, imaginary = 1) * w
  
  #------------------ Strike Exponential -----------------------------------------------------------------------
  k <- (T - S)/100 * (100 - K)
  
  ### Compute Strike Term
  StrikeTerm <- exp(x * k)
  
  #------------------ First Term -----------------------------------------------------------------------
  ### Compute First Term
  FirstTerm <- exp(-(S - t) * rt) * exp(-(T - S) * (1 + x) * rt)
  
  #------------------ Expectation -----------------------------------------------------------------------
  ### Compute Second Expectation
  Expectation <- exp(alpha_qs_ES(x, T - S) + ### initial for alpha
                       alpha_tilde_qs_EtES(x, tau = S - t, T = T, S = S) +
                       beta_tilde_qs_EtES(x, tau = S - t, T = T, S = S) * st)
  
  #------------------ Final Price -----------------------------------------------------------------------
  ### Align dimensions
  StrikeTerm <- as.vector(StrikeTerm)
  #ZCB <- as.vector(ZCB)
  FirstTerm <- as.vector(FirstTerm)
  Expectation <- as.vector(Expectation)
  
  return(StrikeTerm * 1/ZCB * FirstTerm * Expectation)
}


####################################################################################################################
####################################################################################################################
#------------------------------------- Call Price Function ---------------------------------------------------------
####################################################################################################################
####################################################################################################################
price1MSOFROption_Smile <- function(S, t, T, K, ZCB, rt, st, intS_t) { 
  umax <- 50000
  integrand_q_hat <- function(w) {
    Re(q_hat_Smile(w, S = S, t = t, T = T, K = K, ZCB = ZCB, rt = rt, st = st, intS_t = intS_t) / 
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
StrikeSmile <- seq(96-23, 96+23, 1) ### Strikes for the smile
CallSmile <- rep(NA, length(StrikeSmile))
IVSmile <- rep(NA, length(StrikeSmile))
for (i in 1:length(StrikeSmile)) {
  print(i)
  CallSmile[i] <- price1MSOFROption_Smile(S = S, t = 0, T = T, K = StrikeSmile[i], 
                                          ZCB = ZCB_s(t = 0, T = T, rt = r0, st = s0),
                                          rt = r0, st = s0, intS_t = 0)
  IVSmile[i] = IVTest(spot = K_Smile, strike = StrikeSmile[i], t = 0, T = T,
                      ZCB_s(t = 0, T = T, rt = r0, st = s0), CallSmile[i])
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
     ylab = "Black-76 Implicit Volatility", xlab = "Strike", ylim = c(0.07,0.14))
lines(StrikeSmile, IVSmile_MixedExp, col = "steelblue", lwd = 3, lty = 2)
abline(h = 0)
legend("topright", legend = c("Mixed-Exp", "Normal"),
       col = c("steelblue", "black", "#901a1E"), title = "Imp Vol",
       lwd = c(1, 1, 1), lty = c(2, 1, 3), cex=0.7, inset = 0.00)



SmileData_MixedExp <- as.data.frame(cbind(StrikeSmile, CallSmile_MixedExp, IVSmile_MixedExp))
plot(StrikeSmile, IVSmile_MixedExp, type = "p")
abline(h = 0)






plot(StrikeSmile, IVSmile_MixedExp, type = "l", col = "#901a1E", lwd = 3,
     xlab = "Strike", ylab = "Implied Volatility", main = "Black-76 Implied Volatility Across Strikes")
lines(StrikeSmile, IVSmile_Normal, col = "steelblue", lwd = 3)
legend("topright", legend = c("Mixed-Exponential", "Normal"),
       col = c("#901a1E", "steelblue"), title = "Dist. of Unscheduled Jumps",
       lwd = c(3, 3), cex=0.7, inset = 0.00)

plot(StrikeSmile, IVSmile, type = "l", col = "black")
lines(StrikeSmile, IVSmile_MixedExp, col = "#901a1E")
lines(StrikeSmile, IVSmile_Normal, col = "steelblue")

SmileData <- as.data.frame(cbind(StrikeSmile, CallSmile_MixedExp, IVSmile_MixedExp,
                                 CallSmile_Normal, IVSmile_Normal,
                                 CallSmile, IVSmile))
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






