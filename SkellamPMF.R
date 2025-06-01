

n_samples <- 1e6
Skellam <- nu_D * (rpois(n_samples, mu_u * zeta0) - rpois(n_samples, mu_d * zeta0))

### Estimate PMF: round X to exact discrete support points
pmfSkellam <- table(factor(Skellam, levels = sort(unique(Skellam)))) / n_samples
dfSkellam <- data.frame(x = as.numeric(names(pmfSkellam)), prob = as.numeric(pmfSkellam))

### Plot estimated PMF
plot((dfSkellam$x)*1, dfSkellam$prob, type = "h", lwd = 2,
     xlab = "Jump size of target rate (%?)", ylab = "Probability",
     main = "Monte Carlo Estimate of Modified Skellam PMF")
points((dfSkellam$x)*1, dfSkellam$prob, pch = 16, col = "blue")
abline(v = -0.0003, col = "red")


ggplot(dfSkellam, aes(x = x*1)) +
  geom_segment(aes(y = prob, xend = x*1, yend = 0), color = "#666666", linetype = "solid", size = 0.3) +
  geom_point(aes(y = prob), color = "#901a1E", size = 1.3) +
  geom_vline(xintercept = -0.0003, col = "hotpink3", size = 1) +
  labs(title = "Probability Mass Function of Modified Skellam Distribution",
       x = "Jump size of target rate [(%)?]", y = "Probability") +
  #theme(legend.position = "none") +
  theme_minimal()


x_grid <- seq(-0.1, 0.1, by = 0.0001)
plot(x_grid, dnorm(x = x_grid, mean = theta_theta, sd = omega))


#------------------------------

n_samples <- 1e6
XTest <- 1/400 * (rpois(n_samples, 0.6) - rpois(n_samples, 0.1))

# Estimate PMF: round X to exact discrete support points
pmfTest <- table(factor(XTest, levels = sort(unique(XTest)))) / n_samples
dfTest <- data.frame(x = as.numeric(names(pmfTest)), prob = as.numeric(pmfTest))

# Plot estimated PMF
plot((dfTest$x)*100, dfTest$prob, type = "h", lwd = 2,
     xlab = "Jump size of target rate (%)", ylab = "Probability",
     main = "Monte Carlo Estimate of Modified Skellam PMF")
points((dfTest$x)*100, dfTest$prob, pch = 16, col = "blue")



dens <- density(X, bw = "nrd", adjust = 1)

# Plot estimated PDF
plot(dens, main = "Monte Carlo Estimate of PDF of X", xlab = "x", ylab = "Density", col = "darkgreen", lwd = 2)
grid()



####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
ggplot(dfSkellam, aes(x = x*1)) +
  geom_segment(aes(y = prob, xend = x*1, yend = 0), color = "#F0EAD6", linetype = "solid", size = 1) +
  geom_point(aes(y = prob), color = "#4D5D53", size = 1.7) +
  geom_vline(xintercept = -0.0003, col = "#901a1E", size = 1.7, linetype = "dashed") +
  labs(title = "PMF of Modified Skellam Distribution",
       x = "Jump Size of Target Rate", y = "Probability") +
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
  ) +
  #annotate(geom="text", x=0.07, y=0.015, label=paste("Mean"),
  #         color="#901a1E", size = 17) +
  annotate(
    geom = "text",
    x = 0.07, y = 0.0145,
    label = expression(nu^italic(D) * (delta^italic(u) * theta^zeta - delta^italic(d) * theta^zeta)),
    color = "#901a1E", size = 17, parse = TRUE
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
