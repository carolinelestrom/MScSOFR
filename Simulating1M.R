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
jump <- 10/360
jump_ <- 9/360

theta0 <- theta_theta
xi0 <- theta_theta
zeta0 <- theta_zeta
s0 <- 0.003
r0 <- 0.03





####################################################################################################################
####################################################################################################################
#---------------------------------- Calculations -------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Simulation of Master Path
set.seed(NULL) #set.seed(NULL), 23, 333, 3333333
sSim <- rSim <- numeric(nsteps)
sSim[1] <- s0
rSim[1] <- r0
rSim_NN <- rSim_NME <- rSim_SN <- rSim_SME <- numeric(nsteps)
rSim_NN[1] <- rSim_NME[1] <- rSim_SN[1] <- rSim_SME[1] <- r0
J_P_MixedExp <- 0
J_P_Normal <- 0
N_D <- 0
thetam <- theta0
xim <- xi0
J_D_Normal <- 0
zetam <- zeta0
J_D_Skellam <- 0

for (i in 2:nsteps) {
  print(i)
  t <- i * dt
  
  ### Simulate sm at each fine time step
  #sm <- sm + kappa_s_jump * (theta_s_jump - sm) * dt + sqrt(dt) * sigma_s_jump * rnorm(1)
  sSim[i] <- sSim[i - 1] * exp(-kappa_s_jump * dt) +
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
    Z_MixedExp <- ifelse(U < pu_J, Z_u, Z_d)
    Z_Normal <- rnorm(N_P, mean = mu_P, sd = sigma_P)
  } else {
    Z_MixedExp <- 0
    Z_Normal <- 0
  }
  ### Update the jump process
  J_P_MixedExp <- sum(Z_MixedExp) ### Compute compound Poisson sum
  J_P_Normal <- sum(Z_Normal) ### Compute compound Poisson sum
  print(J_P_MixedExp)
  print(J_P_Normal)
  
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
  
  #X_jm <- matrix(data = c(xim, thetam),
  #              nrow = 2, ncol = n???,
  #              byrow = TRUE)
  
  dW_zeta <- rnorm(1)
  zetam <- max(zetam + kappa_zeta * (theta_zeta - zetam) * dt + sigma_zeta * sqrt(pmax(zetam, 0)) * sqrt(dt) * dW_zeta +
                 1 / 4 * sigma_zeta^2 * dt * (dW_zeta^2 - 1), 0)
  dW_xi_d <- rnorm(1)
  #xim_d <- max(xim_d + kappa_d * (theta_d - xim_d) * dt + sigma_d * sqrt(pmax(xim_d, 0)) * sqrt(dt) * dW_xi_d +
  #               1 / 4 * sigma_d^2 * dt * (dW_xi_d^2 - 1), 0)
  #xim_d <- xim_u
  
  ### Generate deterministic jump size
  J_D_Normal <- rnorm(1, mean = gamma_Q + xim, sd = omega) #rnorm(1, mean = gamma_Q + t(Gamma_Q) %*% X_j[, i - 1], sd = omega)
  J_D_Skellam <- nu_D * (rpois(1, lambda = gamma_Q + mu_u * zetam) - rpois(1, lambda = gamma_Q + mu_d * zetam))
  
  if (abs(t - jump) < 1e-5) { #(t == jump) #(time_series[i] %in% jump_dates) #(is_close(t, jump_dates))
    ### Increment the counting process
    N_D <- 1
  } else {
    N_D <- 0
  }
  print(J_D_Normal * N_D)
  print(J_D_Skellam * N_D)
  
  rSim[i] <- rSim[i - 1]  ### Updating at every time step
  rSim_NN[i] <- rSim_NN[i - 1] + J_P_Normal + J_D_Normal * N_D  ### Updating at every time step
  rSim_NME[i] <- rSim_NME[i - 1] + J_P_MixedExp + J_D_Normal * N_D  ### Updating at every time step
  rSim_SN[i] <- rSim_SN[i - 1] + J_P_Normal + J_D_Skellam * N_D  ### Updating at every time step
  rSim_SME[i] <- rSim_SME[i - 1] + J_P_MixedExp + J_D_Skellam * N_D  ### Updating at every time step
}

### Plotting the results
plot(seq(0, 1/12, length.out = nsteps), rSim + sSim, type = "l", ylim = c(0,0.13), lwd = 3)
lines(seq(0, 1/12, length.out = nsteps), rSim_NN + sSim, col = "steelblue", lwd = 3, lty = 2)
lines(seq(0, 1/12, length.out = nsteps), rSim_NME + sSim, col = "#901a1E", lwd = 3, lty = 3)
lines(seq(0, 1/12, length.out = nsteps), rSim_SN + sSim, col = "#39641c", lwd = 3, lty = 2)
lines(seq(0, 1/12, length.out = nsteps), rSim_SME + sSim, col = "#666666", lwd = 3, lty = 3)
abline(v = jump, lwd = 1, lty = 2)
legend("topright",
       legend = expression(
         "No Jumps",
         paste(italic(N), "/", italic(N)),
         paste(italic(N), "/", italic(ME)),
         paste(italic(S), "/", scriptstyle("N")),
         paste(italic(S), "/", scriptstyle("ME")),
         "Jump"
       ),
       col = c("black", "#901a1E", "steelblue", "#666666", "#39641c", "black"),
       title = expression("Dist. of " * J^D * " / " * J^P),
       lwd = c(3, 3, 3, 3, 3, 1), lty = c(1, 1, 1, 1, 1, 2), cex = 0.7, inset = 0.00)



Sim1MDF1 <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "Diffusion1" = rSim + sSim,
                      "NN1" = rSim_NN + sSim, "NME1" = rSim_NME + sSim,
                      "SN1" = rSim_SN + sSim, "SME1" = rSim_SME + sSim)
Sim1MDF2 <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "Diffusion2" = rSim + sSim,
                       "NN2" = rSim_NN + sSim, "NME2" = rSim_NME + sSim,
                       "SN2" = rSim_SN + sSim, "SME2" = rSim_SME + sSim)
Sim1MDF3 <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "Diffusion3" = rSim + sSim,
                       "NN3" = rSim_NN + sSim, "NME3" = rSim_NME + sSim,
                       "SN3" = rSim_SN + sSim, "SME3" = rSim_SME + sSim)
Sim1MDF <- data.frame(Sim1MDF1, Sim1MDF2, Sim1MDF3)
Sim1MDF$time <- Sim1MDF1$t

ggplot(Sim1MDF, aes(x = time)) +
  geom_line(aes(y = NN1, color = "NN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NME1, color = "NME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SN1, color = "SN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SME1, color = "SME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = Diffusion1, color = "Diffusion"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NN2, color = "NN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NME2, color = "NME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SN2, color = "SN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SME2, color = "SME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = Diffusion2, color = "Diffusion"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NN3, color = "NN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NME3, color = "NME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SN3, color = "SN"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = SME3, color = "SME"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = Diffusion3, color = "Diffusion"), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = jump, linetype = "dashed", size = 0.3) +
  labs(title = "Simulations",
       x = "t", y = "") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Distribution", values = c("NN" = "#901a1E", "NME" = "steelblue", 
                                                            "SN" = "#39641c", "SME" = "#666666", 
                                                            "Diffusion" = "hotpink3", "Jump" = "black")) +
  #theme(legend.position = "none") +
  theme_minimal()


ggplot(Sim1MDF, aes(x = time)) +
  geom_line(aes(y = NN1, color = "NN1"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = NME1, color = "NME1"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = SN1, color = "SN1"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = SME1, color = "SME1"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = Diffusion1, color = "Diffusion1"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NN2, color = "NN2"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = NME2, color = "NME2"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = SN2, color = "SN2"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = SME2, color = "SME2"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = Diffusion2, color = "Diffusion2"), linetype = "dashed", size = 0.3) +
  geom_line(aes(y = NN3, color = "NN3"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = NME3, color = "NME3"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = SN3, color = "SN3"), linetype = "solid", size = 0.7) +
  geom_line(aes(y = SME3, color = "SME3"), linetype = "dotted", size = 0.7) +
  geom_line(aes(y = Diffusion3, color = "Diffusion3"), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = jump, linetype = "dashed", size = 0.3) +
  labs(title = "Simulations",
       x = "t", y = "") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Size Distribution",
                     values = c("NN1" = "#901a1E", "NME1" = "lightcoral", 
                                "SN1" = "orchid4", "SME1" = "lightpink2",
                                "Diffusion1" = "hotpink3",
                                "NN2" = "mediumseagreen", "NME2" = "darkolivegreen4", 
                                "SN2" = "#39641c", "SME2" = "forestgreen",
                                "Diffusion2" = "springgreen4",
                                "NN3" = "deepskyblue4", "NME3" = "steelblue", 
                                "SN3" = "cadetblue4", "SME3" = "cyan4",
                                "Diffusion3" = "#0a5963")) +
  #theme(legend.position = "none") +
  theme_minimal()



####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
Sim1MDF <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "Diff" = rSim + sSim,
                       "NN" = rSim_NN + sSim, "NME" = rSim_NME + sSim,
                       "SN" = rSim_SN + sSim, "SME" = rSim_SME + sSim)

ggplot(Sim1MDF, aes(x = t)) +
  geom_line(aes(y = Diff, color = "Diffusion"), linetype = "dashed", size = 1) +
  geom_line(aes(y = NN, color = "N-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = NME, color = "N-ME"), linetype = "dotted", size = 1.7) +
  geom_line(aes(y = SN, color = "S-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = SME, color = "S-ME"), linetype = "dotted", size = 1.7) +
  geom_vline(xintercept = 10/360, col = "#4D5D53", linetype = "dashed", size = 1) +
  labs(title = expression("Simulations of " * r[t]^s * " for " * t %in% "[0, 1/12]"),
       x = "t", y = "") +
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



