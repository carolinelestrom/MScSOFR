

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
         "Jump"
       ),
       col = c("#901a1E", "steelblue", "#666666", "#39641c", "black"),
       title = expression("Dist. of " * J^D * " / " * J^P),
       lwd = c(1, 1, 1, 1, 1), lty = c(1, 1, 1, 1, 2), cex = 0.7, inset = 0.00)
legend("topright",
       legend = c("Normal/Mixed-Exp", "Normal/Normal", "Skellam/Mixed-Exp", "Skellam/Normal", "Scheduled Jump"),
       col = c("#901a1E", "steelblue", "#666666", "#39641c", "black"), title = "Dist. of J_D / J_P",
       lwd = c(1, 1, 1, 1, 1), lty = c(1, 1, 1, 1, 2), cex=0.7, inset = 0.00)



#IVTimeDF_Jump <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "NN" = IV_BS_Normal, "NME" = IV_BS_MixedExp,
#                       "SN" = IV_BS_NormalSkellam, "SME" = IV_BS_MixedExpSkellam,
#                       "Diffusion" = IV_BS_Diffusion)
ggplot(IVTimeDF_Jump, aes(x = t)) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = NN, color = "NN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = NME, color = "NME"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = SN, color = "SN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = SME, color = "SME"), linetype = "dotted", size = 1.3) +
  geom_vline(xintercept = 10/360, col = "black", linetype = "dashed", size = 1) +
  labs(title = "Black-76 Implicit Volatility",
       x = "t", y = "") +
  ylim(0.00, 0.1) +
  scale_color_manual(name = "Jump Distribution", values = c("Diffusion" = "hotpink3", "NN" = "#901a1E", "NME" = "steelblue", "SN" = "#39641c", "SME" = "#666666", "Diffusion" = "hotpink3")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")



IVTimeDF_NoJump <- data.frame("t" = seq(0, 1/12, length.out = nsteps), "NN" = IV_BS_Normal, "NME" = IV_BS_MixedExp,
                            "SN" = IV_BS_NormalSkellam, "SME" = IV_BS_MixedExpSkellam,
                            "Diffusion" = IV_BS_Diffusion)
ggplot(IVTimeDF_NoJump, aes(x = t)) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = NN, color = "NN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = NME, color = "NME"), linetype = "dashed", size = 1.3) +
  geom_line(aes(y = SN, color = "SN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = SME, color = "SME"), linetype = "dashed", size = 1.3) +
  labs(title = "Black-76 Implicit Volatility",
       x = "t", y = "") +
  ylim(0.00, 0.1) +
  scale_color_manual(name = "Jump Distribution", values = c("Diffusion" = "hotpink3", "NN" = "#901a1E", "NME" = "steelblue", "SN" = "#39641c", "SME" = "#666666", "Diffusion" = "hotpink3")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")




####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
### Jump
IVTimeDF_Jump_Cut <- head(IVTimeDF_Jump, 8800) #8400
ggplot(IVTimeDF_Jump_Cut, aes(x = t)) +
  #geom_line(aes(y = Diffusion, color = "Diff"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = NN, color = "N-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = NME, color = "N-ME"), linetype = "dotted", size = 1.7) +
  geom_line(aes(y = SN, color = "S-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = SME, color = "S-ME"), linetype = "dotted", size = 1.7) +
  geom_vline(xintercept = 10/360, col = "#4D5D53", linetype = "dashed", size = 1) +
  labs(title = "Black-76 Implicit Volatility",
       x = "t", y = "IV") +
  ylim(0.00, 0.1) +
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
    legend.background = element_rect(fill = "white", color = "#4D5D53", size = 0.7),
    legend.box.background = element_rect(color = "#4D5D53", size = 0.7)
  )

### No Jump
IVTimeDF_NoJump_Cut <- head(IVTimeDF_NoJump, 8800) #8400
ggplot(IVTimeDF_NoJump_Cut, aes(x = t)) +
  #geom_line(aes(y = Diffusion, color = "Diff"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = NN, color = "N-N"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = NME, color = "N-ME"), linetype = "dotted", size = 1.3) +
  geom_line(aes(y = SN, color = "S-N"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = SME, color = "S-ME"), linetype = "dotted", size = 1.3) +
  labs(title = "Black-76 Implicit Volatility",
       x = "t", y = "IV") +
  ylim(0.00, 0.1) +
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