####################################################################################################################
####################################################################################################################
#--------------------------------- Functions - PDF Optimization ----------------------------------------------------
####################################################################################################################
####################################################################################################################
NormalPDF <- function(x, mean, sd){
  1 / (sqrt(2 * pi * sd^2)) * exp(- ((x - mean)^2) / (2 * sd^2))
}
MixedExpPDF <- function(x, pu_J, qd_J, lambda_J, eta_J, q_J, theta_J){
  if (x >= 0){
    pu_J * sum(lambda_J * eta_J * exp(-eta_J * x))
  } else if (x < 0){
    qd_J * sum(q_J * theta_J * exp(theta_J * x))
  }
}



### Objective function to minimize squared error
SSE <- function(params, x_grid, n, m, mean, sd) {
  ### Extract parameters from vector
  pu_J <- params[1]
  qd_J <- 1 - pu_J  ### Ensures pu_J + qd_J = 1
  
  lambda_J <- params[2:(1 + m)]
  eta_J <- params[(2 + m):(1 + 2 * m)]
  q_J <- params[(2 + 2 * m):(1 + 2 * m + n)]
  theta_J <- params[(2 + 2 * m + n):(1 + 2 * m + 2 * n)]
  
  ### Enforce constraints: sum(lambda_J) = 1, sum(q_J) = 1
  lambda_J <- lambda_J / sum(lambda_J)
  q_J <- q_J / sum(q_J)
  
  ### Enforce positivity of first elements
  lambda_J[1] <- abs(lambda_J[1])
  q_J[1] <- abs(q_J[1])
  
  ### Ensure constraints: sum(lambda_J * eta_J) > 0 and sum(q_J * theta_J) > 0
  if (sum(lambda_J * eta_J) <= 0 || sum(q_J * theta_J) <= 0) {
    return(Inf)
  }
  
  ### Compute normal PDF values
  NormalVals <- sapply(x_grid, function(x) NormalPDF(x, mean, sd))
  
  ### Compute mixed exponential PDF values
  MixedExpVals<- sapply(x_grid, function(x) MixedExpPDF(x, pu_J, qd_J, lambda_J, eta_J, q_J, theta_J))
  
  ### Compute sum of squared errors
  sum((MixedExpVals - NormalVals)^2)
}

### Function to fit the mixed exponential approximation
FitMixedExp <- function(x_grid, n, m, mean, sd) {
  pu_J_init <- 0.5
  ### Initial parameter guesses
  lambda_J_init <- abs(rnorm(m, mean = 0.5, sd = 0.2))
  eta_J_init <- runif(m, 100, 500)
  q_J_init <- abs(rnorm(n, mean = 0.5, sd = 0.2))
  theta_J_init <- runif(n, 100, 500)
  
  init_params <- c(pu_J_init, lambda_J_init, eta_J_init, q_J_init, theta_J_init)
  
  ### Bounds for optimization
  #lower_bounds <- c(0, rep(-Inf, m), rep(100, m), rep(-Inf, n), rep(100, n))
  #upper_bounds <- c(1, rep(Inf, m), rep(1000, m), rep(Inf, n), rep(1000, n))
  lower_bounds <- c(0, rep(0, m), rep(100, m), rep(0, n), rep(100, n))
  upper_bounds <- c(1, rep(1, m), rep(1000, m), rep(1, n), rep(1000, n))
  
  ### Optimize using non-linear least squares
  result <- optim(
    par = init_params,
    fn = SSE,
    x_grid = x_grid,
    n = n,
    m = m,
    mean = mean,
    sd = sd,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds
  )
  
  ### Extract optimized parameters
  params_opt <- result$par
  pu_J_opt <- params_opt[1]
  qd_J_opt <- 1 - pu_J_opt
  lambda_J_opt <- params_opt[2:(1 + m)]
  eta_J_opt <- params_opt[(2 + m):(1 + 2 * m)]
  q_J_opt <- params_opt[(2 + 2 * m):(1 + 2 * m + n)]
  theta_J_opt <- params_opt[(2 + 2 * m + n):(1 + 2 * m + 2 * n)]
  
  ### Normalize lambda_J and q_J
  lambda_J_opt <- lambda_J_opt / sum(lambda_J_opt)
  q_J_opt <- q_J_opt / sum(q_J_opt)
  
  ### Return results
  list(
    pu_J = pu_J_opt,
    qd_J = qd_J_opt,
    lambda_J = lambda_J_opt,
    eta_J = eta_J_opt,
    q_J = q_J_opt,
    theta_J = theta_J_opt,
    sse = result$value
  )
}


n <- 3  ### Number of negative exponentials
m <- 3  ### Number of positive exponentials
x_grid <- seq(-0.1, 0.1, by = 0.0001)
#x_grid <- seq(mu_P - 0.00005, mu_P + 0.00005, length.out = 10000)
x_grid <- seq(mu_P - 0.035, mu_P + 0.035, length.out = 10000)

res <- FitMixedExp(x_grid, n, m, mean = mu_P, sd = (sigma_P + 0.01))

### Print results
print(res)




plot(x_grid, NormalPDF(x_grid, mean = mu_P, sd = sigma_P), type = "l", col = "black", lwd = 3)
lines(x_grid, dnorm(x = x_grid, mean = mu_P, sd = sigma_P), col = "steelblue", lty = 3, lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpPDF(x, pu_J = res$pu_J, qd_J = res$qd_J,
                                                     lambda_J = res$lambda_J, eta_J = res$eta_J,
                                                     q_J = res$q_J, theta_J = res$theta_J)), col = "#901a1E", lty = 3, lwd = 3)


plot(x_grid, NormalPDF(x_grid, mean = mu_P, sd = sigma_P), type = "l", col = "black", lwd = 3)
lines(x_grid, dnorm(x = x_grid, mean = mu_P, sd = sigma_P), col = "steelblue", lty = 3, lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpPDF(x, pu_J = res_CDF$pu_J, qd_J = res_CDF$qd_J,
                                                     lambda_J = res_CDF$lambda_J, eta_J = res_CDF$eta_J,
                                                     q_J = res_CDF$q_J, theta_J = res_CDF$theta_J)), col = "#901a1E", lty = 3, lwd = 3)


####################################################################################################################
####################################################################################################################
#--------------------------------- Functions - CDF Optimization ----------------------------------------------------
####################################################################################################################
####################################################################################################################
NormalCDF <- function(x, mean, sd) {
  pnorm(x, mean, sd)
}

MixedExpCDF <- function(x, pu_J, qd_J, lambda_J, eta_J, q_J, theta_J) {
  if (x >= 0) {
    return(qd_J * sum(q_J * exp(theta_J * 0)) +
             pu_J * sum(lambda_J * (1 - exp(-eta_J * x))))
  } else if (x < 0) {
    return(qd_J * sum(q_J * exp(theta_J * x)))
  }
}

### Objective function to minimize squared error in CDFs
SSE_CDF <- function(params, x_grid, n, m, mean, sd) {
  pu_J <- params[1]
  qd_J <- 1 - pu_J  ### Ensures pu_J + qd_J = 1
  
  lambda_J <- params[2:(1 + m)]
  eta_J <- params[(2 + m):(1 + 2 * m)]
  q_J <- params[(2 + 2 * m):(1 + 2 * m + n)]
  theta_J <- params[(2 + 2 * m + n):(1 + 2 * m + 2 * n)]
  
  ### Enforce constraints: sum(lambda_J) = 1, sum(q_J) = 1
  lambda_J <- lambda_J / sum(lambda_J)
  q_J <- q_J / sum(q_J)
  
  ### Ensure positivity
  lambda_J[1] <- abs(lambda_J[1])
  q_J[1] <- abs(q_J[1])
  
  ### Compute normal CDF values
  NormalVals <- sapply(x_grid, function(x) NormalCDF(x, mean, sd))
  
  ### Compute mixed exponential CDF values
  MixedExpVals <- sapply(x_grid, function(x) MixedExpCDF(x, pu_J, qd_J, lambda_J, eta_J, q_J, theta_J))
  
  ### Compute sum of squared errors
  sum((MixedExpVals - NormalVals)^2)
}

### Function to fit the mixed exponential approximation using CDF
FitMixedExp_CDF <- function(x_grid, n, m, mean, sd) {
  pu_J_init <- 0.5
  lambda_J_init <- abs(rnorm(m, mean = 0.5, sd = 0.2))
  eta_J_init <- runif(m, 100, 1000)
  q_J_init <- abs(rnorm(n, mean = 0.5, sd = 0.2))
  theta_J_init <- runif(n, 100, 1000)
  
  init_params <- c(pu_J_init, lambda_J_init, eta_J_init, q_J_init, theta_J_init)
  
  lower_bounds <- c(0, rep(0.001, m), rep(1.001, m), rep(0.001, n), rep(0.001, n))
  upper_bounds <- c(1, rep(10, m), rep(1000, m), rep(10, n), rep(1000, n))
  
  ### Optimize using CDF
  result <- optim(
    par = init_params,
    fn = SSE_CDF,
    x_grid = x_grid,
    n = n,
    m = m,
    mean = mean,
    sd = sd,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds
  )
  
  ### Extract optimized parameters
  params_opt <- result$par
  pu_J_opt <- params_opt[1]
  qd_J_opt <- 1 - pu_J_opt
  lambda_J_opt <- params_opt[2:(1 + m)]
  eta_J_opt <- params_opt[(2 + m):(1 + 2 * m)]
  q_J_opt <- params_opt[(2 + 2 * m):(1 + 2 * m + n)]
  theta_J_opt <- params_opt[(2 + 2 * m + n):(1 + 2 * m + 2 * n)]
  
  ### Normalize lambda_J and q_J
  lambda_J_opt <- lambda_J_opt / sum(lambda_J_opt)
  q_J_opt <- q_J_opt / sum(q_J_opt)
  
  ### Return results
  list(
    pu_J = pu_J_opt,
    qd_J = qd_J_opt,
    lambda_J = lambda_J_opt,
    eta_J = eta_J_opt,
    q_J = q_J_opt,
    theta_J = theta_J_opt,
    sse = result$value
  )
}

### Run the CDF-based approximation
n <- 3  ### Number of negative exponentials
m <- 3  ### Number of positive exponentials
#x_grid <- seq(mu_P - 0.00005, mu_P + 0.00005, length.out = 10000)
x_grid <- seq(mu_P - 0.035, mu_P + 0.035, length.out = 10000)

res_CDF <- FitMixedExp_CDF(x_grid, n, m, mean = mu_P, sd = sigma_P)


### Print results
print(res_CDF)




plot(x_grid, NormalCDF(x_grid, mean = mu_P, sd = sigma_P), type = "l", col = "black", lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpCDF(x, pu_J = res_CDF$pu_J, qd_J = res_CDF$qd_J,
                                                     lambda_J = res_CDF$lambda_J, eta_J = res_CDF$eta_J,
                                                     q_J = res_CDF$q_J, theta_J = res_CDF$theta_J)),
      col = "#901a1E", lty = 3, lwd = 3)

plot(x_grid, NormalCDF(x_grid, mean = mu_P, sd = (sigma_P)), type = "l", col = "black", lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpCDF(x, pu_J = res$pu_J, qd_J = res$qd_J,
                                                     lambda_J = res$lambda_J, eta_J = res$eta_J,
                                                     q_J = res$q_J, theta_J = res$theta_J)),
      col = "#901a1E", lty = 3, lwd = 3)


MEApproxN <- data.frame("x" = x_grid,
                        "N" = NormalCDF(x_grid, mean = mu_P, sd = sigma_P),
                        "ME" = sapply(x_grid, function(x) MixedExpCDF(x, pu_J = res_CDF$pu_J, qd_J = res_CDF$qd_J,
                                                                      lambda_J = res_CDF$lambda_J, eta_J = res_CDF$eta_J,
                                                                      q_J = res_CDF$q_J, theta_J = res_CDF$theta_J)))
ggplot(MEApproxN, aes(x = x)) +
  geom_line(aes(y = N, color = "N"), linetype = "solid", size = 3) +
  geom_line(aes(y = ME, color = "ME"), linetype = "dotted", size = 3) +
  labs(title = expression("Approximate the CDF of" ~ N(mu^P, sigma^P)),
       x = "x", y = "CDF") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Distribution", values = c("N" = "#666666", "ME" = "#901a1E")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")


#--------------------------------------------------------------------

### Run the CDF-based approximation
n <- 3  ### Number of negative exponentials
m <- 3  ### Number of positive exponentials
x_grid <- seq(-3, 3, length.out = 10000)
resTest <- FitMixedExp_CDF(x_grid, n, m, mean = mu_P, sd = 1)
print(resTest)

plot(x_grid, NormalCDF(x_grid, mean = 0, sd = 1), type = "l", col = "black", lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpCDF(x, pu_J = resTest$pu_J, qd_J = resTest$qd_J,
                                                     lambda_J = resTest$lambda_J, eta_J = resTest$eta_J,
                                                     q_J = resTest$q_J, theta_J = resTest$theta_J)),
      col = "#901a1E", lty = 3, lwd = 3)


x_grid <- seq(-0.035, 0.035, length.out = 10000)
plot(x_grid, NormalCDF(x_grid, mean = 0, sd = 0.01), type = "l", col = "black", lwd = 3)
lines(x_grid, sapply(x_grid, function(x) MixedExpCDF(x, pu_J = 0.5, qd_J = 0.5,
                                                     lambda_J = c(8.7303, 2.1666, -10), eta_J = c(213.0215, 236.0406, 237.1139),
                                                     q_J = c(0.0622, 0.0409), theta_J = c(939.7441, 939.8021))),
      col = "#901a1E", lty = 3, lwd = 3)




####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
ggplot(MEApproxN, aes(x = x)) +
  geom_line(aes(y = N, color = "N"), linetype = "solid", size = 3) +
  geom_line(aes(y = ME, color = "ME"), linetype = "dotted", size = 3) +
  labs(title = expression("Approximate CDF of " ~ N(mu^P, sigma^P) ~ " by " ~ ME(eta^P, theta^P)),
       x = "x", y = "CDF") +
  #ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Size Distribution", values = c("N" = "#4D5D53", "ME" = "#901a1E")) +
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
    legend.background = element_rect(fill = "white", color = "#4D5D53", size = 0.7),
    legend.box.background = element_rect(color = "#4D5D53", size = 0.7)
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
"#901a1E", "#39641c", "#666666", "#ffbd38", "#0a5963", "#122947", "#425570"
### Nice
"hotpink3", "orchid4", "steelblue"



