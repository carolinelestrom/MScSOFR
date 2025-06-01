####################################################################################################################
####################################################################################################################
#------------------------------------ Settings ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Parameters
T <- 12/12
S <- 0
nsteps <- 1000
dt <- T / nsteps
Daily <- seq(1, T * 360) / 360  ### Daily steps

### Initial Values
theta0 <- theta_theta
xi0 <- theta_theta
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
#------------------------------------ SOFR ZCB Price ---------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### SOFR Spread ----------------------------------------------------------------------------------------------------
####################################################################################################################
### Exact solution for beta_s
beta_s <- function(tau){
  1/kappa_s_jump * (exp(-kappa_s_jump * tau) - 1)
}

### Exact solution for alpha_s
alpha_s <- function(tau){
  -theta_s_jump * tau + (sigma_s_jump^2 * tau)/(2 * kappa_s_jump^2) +
    (theta_s_jump/kappa_s_jump - (sigma_s_jump^2)/(kappa_s_jump^3)) * (1 - exp(-kappa_s_jump * tau)) +
    (sigma_s_jump^2 * (1 - exp(-2 * kappa_s_jump * tau)))/(4 * kappa_s_jump^3)
}

ZCB_s <- function(t, T, rt, st){
  ### Compute ZCB price
  ZCB <- exp(-(T - t) * rt +  
               alpha_s(T - t) + 
               beta_s(T - t) * st)
  
  return(ZCB)
}

fYTM_s <- function(t, T, rt, st) {
  -log(ZCB_s(t, T, rt, st)) / (T - t)
}




### Normal Scheduled Jumps and Normal Unscheduled Jumps ------------------------------------------------------------
####################################################################################################################
### Zero-Coupon Bond Price Function
ZCB_s_NN <- function(t, T, thetat, xit, rt, st){
  ### Scheduled Jump Contribution
  NormalScheduledJump <- exp(-(T - jumps[1]) * MeanXi(t, jumps_[1], thetat, xit)) *                           ### 1st jump
    exp((-(T - jumps[1]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(t, jumps_[1]) %*% Gamma_Q)) *                ### 1st jump
    exp(-(T - jumps[2]) * MeanXi(t, jumps_[2], thetat, xit)) *                                                ### 2nd jump
    exp((-(T - jumps[2]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[1], jumps_[2]) %*% Gamma_Q)) *        ### 2nd jump
    exp(-(T - jumps[3]) * MeanXi(t, jumps_[3], thetat, xit)) *                                                ### 3rd jump
    exp((-(T - jumps[3]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[2], jumps_[3]) %*% Gamma_Q)) *        ### 3rd jump
    exp(-(T - jumps[4]) * MeanXi(t, jumps_[4], thetat, xit)) *                                                ### 4th jump
    exp((-(T - jumps[4]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[3], jumps_[4]) %*% Gamma_Q)) *        ### 4th jump
    exp(-(T - jumps[5]) * MeanXi(t, jumps_[5], thetat, xit)) *                                                ### 5th jump
    exp((-(T - jumps[5]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[4], jumps_[5]) %*% Gamma_Q)) *        ### 5th jump
    exp(-(T - jumps[6]) * MeanXi(t, jumps_[6], thetat, xit)) *                                                ### 6th jump
    exp((-(T - jumps[6]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[5], jumps_[6]) %*% Gamma_Q)) *        ### 6th jump
    exp(-(T - jumps[7]) * MeanXi(t, jumps_[7], thetat, xit)) *                                                ### 7th jump
    exp((-(T - jumps[7]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[6], jumps_[7]) %*% Gamma_Q)) *        ### 7th jump
    exp(-(T - jumps[8]) * MeanXi(t, jumps_[8], thetat, xit)) *                                                ### 8th jump in reference period
    exp((-(T - jumps[8]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[7], jumps_[8]) %*% Gamma_Q))          ### 8th jump in reference period
  
  ### Define func_qJ_ES to compute alpha_qJ using the ODE
  func_N <- function(tau, alpha_N) {
    nu_P * (exp(mu_P * (-tau) + (sigma_P * (-tau))^2 / 2) - 1)
  }
  
  ### Use a simple RK4 or your own integration routine (e.g., from pracma or deSolve)
  taus <- seq(0, T - t, length.out = 100)
  alpha_N_vals <- numeric(length(taus))
  alpha_N_vals[1] <- 0                       ### Initial value
  dtau <- taus[2] - taus[1]                  ### Time step difference
  
  for (i in 2:length(taus)) {
    #k1 <- dtau * func_N(taus[i - 1], alpha_N_vals[i - 1])
    #k2 <- dtau * func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + k1 / 2)
    #k3 <- dtau * func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + k2 / 2)
    #k4 <- dtau * func_N(taus[i - 1] + dtau, alpha_N_vals[i - 1] + k3)
    #alpha_N_vals[i] <- alpha_N_vals[i - 1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    k1 <- func_N(taus[i - 1], alpha_N_vals[i - 1])
    k2 <- func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + dtau * k1 / 2)
    k3 <- func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + dtau * k2 / 2)
    k4 <- func_N(taus[i - 1] + dtau, alpha_N_vals[i - 1] + dtau * k3)
    alpha_N_vals[i] <- alpha_N_vals[i - 1] + (dtau / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }
  
  alpha_N_val <- tail(alpha_N_vals, 1)             ### Unscheduled Jump Contribution
  
  ### Compute ZCB price
  ZCB <- exp(-(T - t) * rt +  
               alpha_N_val +
               alpha_s(T - t) + 
               beta_s(T - t) * st) *
    NormalScheduledJump
  
  return(ZCB)
}

fYTM_s_NN <- function(t, T, thetat, xit, rt, st) {
  -log(ZCB_s_NN(t, T, thetat, xit, rt, st)) / (T - t)
}



### Skellam Scheduled Jumps and Normal Unscheduled Jumps -----------------------------------------------------------
####################################################################################################################
### Zero-Coupon Bond Price Function
ZCB_s_SN <- function(t, T, zetat, rt, st){
  ### Scheduled Jump Contribution
  SkellamScheduledJump <- exp(alpha_tau8(tau = jumps_[8] - jumps_[7]) + 
                                alpha_tau7(tau = jumps_[7] - jumps_[6]) + 
                                alpha_tau6(tau = jumps_[6] - jumps_[5]) + 
                                alpha_tau5(tau = jumps_[5] - jumps_[4]) + 
                                alpha_tau4(tau = jumps_[4] - jumps_[3]) + 
                                alpha_tau3(tau = jumps_[3] - jumps_[2]) + 
                                alpha_tau2(tau = jumps_[2] - jumps_[1]) + 
                                alpha_tau1(tau = jumps_[1] - t) + 
                                beta_tau1(tau = jumps_[1] - t) * zetat)
  
  ### Define func_qJ_ES to compute alpha_qJ using the ODE
  func_N <- function(tau, alpha_N) {
    nu_P * (exp(mu_P * (-tau) + (sigma_P * (-tau))^2 / 2) - 1)
  }
  
  ### Use a simple RK4 or your own integration routine (e.g., from pracma or deSolve)
  taus <- seq(0, T - t, length.out = 100)
  alpha_N_vals <- numeric(length(taus))
  alpha_N_vals[1] <- 0                       ### Initial value
  dtau <- taus[2] - taus[1]                  ### Time step difference
  
  for (i in 2:length(taus)) {
    #k1 <- dtau * func_N(taus[i - 1], alpha_N_vals[i - 1])
    #k2 <- dtau * func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + k1 / 2)
    #k3 <- dtau * func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + k2 / 2)
    #k4 <- dtau * func_N(taus[i - 1] + dtau, alpha_N_vals[i - 1] + k3)
    #alpha_N_vals[i] <- alpha_N_vals[i - 1] + (k1 + 2*k2 + 2*k3 + k4) / 6
    k1 <- func_N(taus[i - 1], alpha_N_vals[i - 1])
    k2 <- func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + dtau * k1 / 2)
    k3 <- func_N(taus[i - 1] + dtau/2, alpha_N_vals[i - 1] + dtau * k2 / 2)
    k4 <- func_N(taus[i - 1] + dtau, alpha_N_vals[i - 1] + dtau * k3)
    alpha_N_vals[i] <- alpha_N_vals[i - 1] + (dtau / 6) * (k1 + 2*k2 + 2*k3 + k4)
  }
  
  alpha_N_val <- tail(alpha_N_vals, 1)             ### Unscheduled Jump Contribution
  
  ### Compute ZCB price
  ZCB <- exp(-(T - t) * rt +  
               alpha_N_val +
               alpha_s(T - t) + 
               beta_s(T - t) * st) *
    SkellamScheduledJump
  
  return(ZCB)
}

fYTM_s_SN <- function(t, T, zetat, rt, st) {
  -log(ZCB_s_SN(t, T, zetat, rt, st)) / (T - t)
}


### Normal Scheduled Jumps and Mixed-Exponential Unscheduled Jumps -------------------------------------------------
####################################################################################################################
### Exact solution for alpha^E
alpha_E <- function(tau){
  Pos <- numeric(length(lambda_J))
  for(i in 1:length(lambda_J)){
    Pos[i] <- pu_J * lambda_J[i] * eta_J[i] *
      log((tau + eta_J[i]) / (eta_J[i]))
  }
  Neg <- numeric(length(q_J))
  for(i in 1:length(q_J)){
    Neg[i] <- qd_J * q_J[i] * theta_J[i] *
      log((-tau + theta_J[i]) / (theta_J[i]))
  }
  return(nu_P * (-tau + sum(Pos) - sum(Neg)))
}

### Zero-Coupon Bond Price Function
ZCB_s_NME <- function(t, T, thetat, xit, rt, st){
  ### Scheduled Jump Contribution
  NormalScheduledJump <- exp(-(T - jumps[1]) * MeanXi(t, jumps_[1], thetat, xit)) *                           ### 1st jump
    exp((-(T - jumps[1]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(t, jumps_[1]) %*% Gamma_Q)) *                ### 1st jump
    exp(-(T - jumps[2]) * MeanXi(t, jumps_[2], thetat, xit)) *                                                ### 2nd jump
    exp((-(T - jumps[2]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[1], jumps_[2]) %*% Gamma_Q)) *        ### 2nd jump
    exp(-(T - jumps[3]) * MeanXi(t, jumps_[3], thetat, xit)) *                                                ### 3rd jump
    exp((-(T - jumps[3]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[2], jumps_[3]) %*% Gamma_Q)) *        ### 3rd jump
    exp(-(T - jumps[4]) * MeanXi(t, jumps_[4], thetat, xit)) *                                                ### 4th jump
    exp((-(T - jumps[4]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[3], jumps_[4]) %*% Gamma_Q)) *        ### 4th jump
    exp(-(T - jumps[5]) * MeanXi(t, jumps_[5], thetat, xit)) *                                                ### 5th jump
    exp((-(T - jumps[5]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[4], jumps_[5]) %*% Gamma_Q)) *        ### 5th jump
    exp(-(T - jumps[6]) * MeanXi(t, jumps_[6], thetat, xit)) *                                                ### 6th jump
    exp((-(T - jumps[6]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[5], jumps_[6]) %*% Gamma_Q)) *        ### 6th jump
    exp(-(T - jumps[7]) * MeanXi(t, jumps_[7], thetat, xit)) *                                                ### 7th jump
    exp((-(T - jumps[7]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[6], jumps_[7]) %*% Gamma_Q)) *        ### 7th jump
    exp(-(T - jumps[8]) * MeanXi(t, jumps_[8], thetat, xit)) *                                                ### 8th jump in reference period
    exp((-(T - jumps[8]))^2 / 2 * (omega^2 + t(Gamma_Q) %*% VarX(jumps_[7], jumps_[8]) %*% Gamma_Q))          ### 8th jump in reference period
  
  ### Compute ZCB price
  ZCB <- exp(-(T - t) * rt + 
               alpha_E(T - t) + 
               alpha_s(T - t) + 
               beta_s(T - t) * st) *
    NormalScheduledJump
  
  return(ZCB)
}

fYTM_s_NME <- function(t, T, thetat, xit, rt, st) {
  -log(ZCB_s_NME(t, T, thetat, xit, rt, st)) / (T - t)
}

### Skellam Scheduled Jumps and Mixed-Exponential Unscheduled Jumps ------------------------------------------------
####################################################################################################################
### Zero-Coupon Bond Price Function
ZCB_s_SME <- function(t, T, zetat, rt, st){
  ### Scheduled Jump Contribution
  SkellamScheduledJump <- exp(alpha_tau8(tau = jumps_[8] - jumps_[7]) + 
                                alpha_tau7(tau = jumps_[7] - jumps_[6]) + 
                                alpha_tau6(tau = jumps_[6] - jumps_[5]) + 
                                alpha_tau5(tau = jumps_[5] - jumps_[4]) + 
                                alpha_tau4(tau = jumps_[4] - jumps_[3]) + 
                                alpha_tau3(tau = jumps_[3] - jumps_[2]) + 
                                alpha_tau2(tau = jumps_[2] - jumps_[1]) + 
                                alpha_tau1(tau = jumps_[1] - t) + 
                                beta_tau1(tau = jumps_[1] - t) * zetat)
  
  ### Compute ZCB price
  ZCB <- exp(-(T - t) * rt + 
               alpha_E(T - t) + 
               alpha_s(T - t) + 
               beta_s(T - t) * st) *
    SkellamScheduledJump
  
  return(ZCB)
}

fYTM_s_SME <- function(t, T, zetat, rt, st) {
  -log(ZCB_s_SME(t, T, zetat, rt, st)) / (T - t)
}

####################################################################################################################
####################################################################################################################
#----------------------------------- Calculations ------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### ZCB Prices
ZCB_T <- ZCB_T_NN <- ZCB_T_NME <- ZCB_T_SN <- ZCB_T_SME <- numeric(nsteps)
for (i in 1:nsteps) {
  print(i)
  Mat <- i * dt
  ZCB_T[i] <- ZCB_s(0, Mat, r0, s0)
  ZCB_T_NN[i] <- ZCB_s_NN(0, Mat, theta0, xi0, r0, s0)
  ZCB_T_NME[i] <- ZCB_s_NME(0, Mat, theta0, xi0, r0, s0)
  ZCB_T_SN[i] <- ZCB_s_SN(0, Mat, zeta0, r0, s0)
  ZCB_T_SME[i] <- ZCB_s_SME(0, Mat, zeta0, r0, s0)
}


T <- 12/12
nsteps <- 1000
plot(seq(0, T, length.out = nsteps), ZCB_T, type = "l", lwd = 3,
     xlab = "T", ylab = "ZCB", main = "Bond Prices for T in [0, 12/12]")
lines(seq(0, T, length.out = nsteps), ZCB_T_NN, col = "steelblue", lwd = 3, lty = 2)
lines(seq(0, T, length.out = nsteps), ZCB_T_NME, col = "#901a1E", lwd = 3, lty = 3)
lines(seq(0, T, length.out = nsteps), ZCB_T_SN, col = "#39641c", lwd = 3, lty = 2)
lines(seq(0, T, length.out = nsteps), ZCB_T_SME, col = "#666666", lwd = 3, lty = 3)
legend("topright",
       legend = expression(
         "Diffusion",
         paste(italic(N), "/", italic(N)),
         paste(italic(N), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         paste(italic(S), "/", scriptstyle("ME"))
       ),
       col = c("black", "steelblue", "#901a1E", "#39641c", "#666666"),
       title = expression("Dist. of " * J^D * " / " * J^P),
       lwd = c(3, 3, 3, 3, 3), lty = c(1, 2, 3, 2, 3), cex = 0.7, inset = 0.00)



ZCBDF <- data.frame("T" = seq(0, T, length.out = nsteps), "Diffusion" = ZCB_T,
                    "NN" = ZCB_T_NN, "NME" = ZCB_T_NME,
                    "SN" = ZCB_T_SN, "SME" = ZCB_T_SME)

ggplot(ZCBDF, aes(x = T)) +
  geom_line(aes(y = NN, color = "NN"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = NME, color = "NME"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = SN, color = "SN"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = SME, color = "SME"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dashed", size = 0.7) +
  labs(title = "Bond Prices for T in [0, 12/12]",
       x = "T", y = "") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Distribution", values = c("NN" = "#901a1E", "NME" = "steelblue", 
                                                            "SN" = "#39641c", "SME" = "#666666", 
                                                            "Diffusion" = "hotpink3", "Jump" = "black")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")

####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
ggplot(ZCBDF, aes(x = T)) +
  geom_line(aes(y = NN, color = "N-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = NME, color = "N-ME"), linetype = "dotted", size = 1.7) +
  geom_line(aes(y = SN, color = "S-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = SME, color = "S-ME"), linetype = "dotted", size = 1.7) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dashed", size = 1.7) +
  labs(title = expression("Bond Prices for " * t %in% "[0, 12/12]"),
       x = "T", y = "") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Size Distribution", values = c("Diff" = "#4D5D53",
                                                                 "N-N" = "#901a1E", 
                                                                 "N-ME" = "#FFBCD9", 
                                                                 "S-N" = "#39641c", 
                                                                 "S-ME" = "#FF8C00")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(
    plot.title = element_text(size = 47, hjust = 0.5, vjust = 1.5),
    axis.title = element_text(size = 37),
    axis.text.x = element_text(size = 23),
    axis.text.y = element_text(size = 23),
    legend.position = "top",
    legend.text = element_text(size = 23),
    legend.title = element_text(size = 27),
    legend.background = element_rect(fill = "white", color = "black", size = 0.7),
    legend.box.background = element_rect(color = "black", size = 0.7)
  )




####################################################################################################################
####################################################################################################################
#----------------------------------- Color Pallette ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Black/White ish
"#4D5D53", "#F0EAD6"
### Bielefeld Vibes
"#FF8C00", "#404080", "#87A96B", "#FFBCD9", "#996666"
### KU
#"#901a1E", "#39641c", "#666666", "#ffbd38", "#0a5963", "#122947", "#425570"
### Nice
#"hotpink3", "orchid4", "steelblue"



####################################################################################################################
####################################################################################################################
#-------------------------------------- Yields ---------------------------------------------------------------------
####################################################################################################################
####################################################################################################################

plot(seq(0, T, length.out = nsteps)[1:200], ZCB_T[1:200], type = "l", lwd = 3,
     xlab = "T", ylab = "ZCB", main = "Bond Prices for T in [0, 12/12]")
lines(seq(0, T, length.out = nsteps)[1:200], ZCB_T_NN[1:200], col = "steelblue", lwd = 3, lty = 2)
lines(seq(0, T, length.out = nsteps)[1:200], ZCB_T_NME[1:200], col = "#901a1E", lwd = 3, lty = 3)
lines(seq(0, T, length.out = nsteps)[1:200], ZCB_T_SN[1:200], col = "#39641c", lwd = 3, lty = 2)
lines(seq(0, T, length.out = nsteps)[1:200], ZCB_T_SME[1:200], col = "#666666", lwd = 3, lty = 3)
legend("topright",
       legend = expression(
         "Diffusion",
         paste(italic(N), "/", italic(N)),
         paste(italic(N), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         paste(italic(S), "/", scriptstyle("ME"))
       ),
       col = c("black", "steelblue", "#901a1E", "#39641c", "#666666"),
       title = expression("Dist. of " * J^D * " / " * J^P),
       lwd = c(3, 3, 3, 3, 3), lty = c(1, 2, 3, 2, 3), cex = 0.7, inset = 0.00)


### Yields
Yield <- Yield_N <- Yield_ME <- numeric(nsteps)
for (i in 1:nsteps) {
  print(i)
  Mat <- i * dt
  Yield[i] <- fYTM_s(0, Mat, r0, s0)
  Yield_N[i] <- fYTM_s_N(0, Mat, r0, s0)
  Yield_ME[i] <- fYTM_s_ME(0, Mat, r0, s0)
}
plot(seq(0, T, length.out = nsteps), Yield)
lines(seq(0, T, length.out = nsteps), Yield_N, col = "steelblue")
lines(seq(0, T, length.out = nsteps), Yield_ME, col = "#901a1E")






