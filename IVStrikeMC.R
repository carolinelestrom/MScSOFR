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
paths <- 10000  ### Number of simulated paths

### Initial Values
theta0 <- 0
xi0 <- 0
s0 <- 0.00 #0.003
r0 <- 0.03 #0.03
st <- rep(r0 + s0, paths)
int_tT <- rep(0, paths)
int_ST <- rep(0, paths)


####################################################################################################################
####################################################################################################################
#------------------------------------ Simulations ------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Vasicek Process Simulation (Exact or Euler)
set.seed(NULL)
for (i in 1:nsteps) {
  print(i)
  t <- i * dt
  
  #sm <- sm + kappa_s_jump * (theta_s_jump - sm) * dt + sqrt(dt) * sigma_s_jump * rnorm(paths)
  st <- st * exp(-kappa_s_jump * dt) +
    theta_s_jump * (1 - exp(-kappa_s_jump * dt)) +
    sigma_s_jump * sqrt((1 - exp(-2 * kappa_s_jump * dt)) / (2 * kappa_s_jump)) * rnorm(paths)
  
  int_tT <- int_tT + st * dt
  
  if (t >= S) {
    int_ST <- int_ST + st * dt
  }
}


####################################################################################################################
####################################################################################################################
#-------------------------------- Monte Carlo Option Pricing -------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Strikes
StrikeSmile <- seq(96-23, 96+23, 1) ### Strikes for the smile
kSmile <- (T - S) * (100 - StrikeSmile) / 100

### Monte Carlo Option Pricing
MCPrice <- numeric(length(StrikeSmile))
for (i in 1:length(StrikeSmile)){
  MCPrice[i] <- (100 * ZCB_s_Normal(t = 0, T = T, rt = r0, st = s0) / (T - S)) *
    mean(exp(-int_tT) * pmax(kSmile[i] - int_ST, 0))
}


####################################################################################################################
####################################################################################################################
#-------------------------------- Implicit Volatility --------------------------------------------------------------
####################################################################################################################
####################################################################################################################
IVSmile_MC <- numeric(length(StrikeSmile))
for (i in 1:length(StrikeSmile)) {
  print(i)
  IVSmile_MC[i] = IVTest(spot = K_Smile_Normal, strike = StrikeSmile[i],
                         t = 0, T = T, ZCB_s_Normal(t = 0, T = T, rt = r0, st = s0), MCPrice[i])
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

