####################################################################################################################
####################################################################################################################
#------------------------------------- Moments Functions -----------------------------------------------------------
####################################################################################################################
####################################################################################################################
MeanXi <- function(t, T, thetat, xit){
  exp(-kappa_xi * (T - t)) * xit + 
    theta_theta * (1 - exp(-kappa_xi * (T - t))) +
    (thetat - theta_theta) * kappa_xi/(kappa_xi - kappa_theta) * 
    (exp(-kappa_theta * (T - t)) - exp(-kappa_xi * (T - t)))
}

MeanTheta <- function(t, T, thetat){
  exp(-kappa_theta * (T - t)) * thetat + theta_theta * (1 - exp(-kappa_theta * (T - t)))
}

MeanX <- function(t, T, thetat, xit){
  matrix(data = c(exp(-kappa_xi * (T - t)) * xit + 
                    theta_theta * (1 - exp(-kappa_xi * (T - t))) +
                    (thetat - theta_theta) * kappa_xi/(kappa_xi - kappa_theta) * 
                    (exp(-kappa_theta * (T - t)) - exp(-kappa_xi * (T - t))),
                  exp(-kappa_theta * (T - t)) * thetat + 
                    (1 - exp(-kappa_theta * (T - t))) * theta_theta),
         byrow = TRUE, nrow = 2)
}

VarX <- function(t, T){
  matrix(data = c((sigma_xi^2 / (2 * kappa_xi) - 
                     2 * rho * sigma_xi * sigma_theta * kappa_xi / (kappa_xi - kappa_theta) * 1 / (2 * kappa_xi) +
                     sigma_theta^2 * kappa_xi^2 / ((kappa_xi - kappa_theta)^2) * 1 / (2 * kappa_xi)) * 
                    (1 - exp(-2 * kappa_xi * (T - t))) +
                    sigma_theta^2 * kappa_xi^2 / ((kappa_xi - kappa_theta)^2) * 1 / (2 * kappa_theta) *
                    (1 - exp(-2 * kappa_theta * (T - t))) +
                    (2 * rho * sigma_xi * sigma_theta * kappa_xi / (kappa_xi - kappa_theta) * 1 / (kappa_xi + kappa_theta) -
                       2 * sigma_theta^2 * kappa_xi^2 / ((kappa_xi - kappa_theta)^2) * 1 / (kappa_xi + kappa_theta)) *
                    (1 - exp(-(kappa_xi + kappa_theta) * (T - t))),
                  ((rho * sigma_xi * sigma_theta) / (kappa_xi + kappa_theta) - 
                     sigma_theta^2 * kappa_xi / (kappa_xi - kappa_theta) * 1 / (kappa_xi + kappa_theta)) *
                    (1 - exp(-(kappa_xi + kappa_theta) * (T - t))),
                  ((rho * sigma_xi * sigma_theta) / (kappa_xi + kappa_theta) - 
                     sigma_theta^2 * kappa_xi / (kappa_xi - kappa_theta) * 1 / (kappa_xi + kappa_theta)) *
                    (1 - exp(-(kappa_xi + kappa_theta) * (T - t))), 
                  sigma_theta^2/(2 * kappa_theta) * (1 - exp(-2 * kappa_theta * (T - t)))),
         byrow = TRUE, nrow = 2)
}






