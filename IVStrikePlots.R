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




IVSmileDF <- data.frame("K" = StrikeSmile, "NN" = IVSmile_Normal, "NME" = IVSmile_MixedExp,
                        "SN" = IVSmile_NormalSkellam, "SME" = IVSmile_MixedExpSkellam, "Diffusion" = IVSmile)
ggplot(IVSmileDF, aes(x = StrikeSmile)) +
  geom_line(aes(y = NN, color = "NN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = NME, color = "NME"), linetype = "dashed", size = 1.3) +
  geom_line(aes(y = SN, color = "SN"), linetype = "solid", size = 1.3) +
  geom_line(aes(y = SME, color = "SME"), linetype = "dashed", size = 1.3) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dotted", size = 1.3) +
  labs(title = "Black-76 Implicit Volatility",
       x = "Strike", y = "IV") +
  ylim(0.05, 0.25) +
  scale_color_manual(name = "Jump Distribution", values = c("NN" = "#901a1E", "NME" = "steelblue", "SN" = "#39641c", "SME" = "#666666", "Diffusion" = "hotpink3")) +
  #theme(legend.position = "none") +
  theme_minimal() %+replace%
  theme(legend.position = "top")


####################################################################################################################
####################################################################################################################
#--------------------------------------- Final Plot ----------------------------------------------------------------
####################################################################################################################
####################################################################################################################
ggplot(IVSmileDF, aes(x = StrikeSmile)) +
  geom_line(aes(y = NN, color = "N-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = NME, color = "N-ME"), linetype = "dashed", size = 1.7) +
  geom_line(aes(y = SN, color = "S-N"), linetype = "solid", size = 1.7) +
  geom_line(aes(y = SME, color = "S-ME"), linetype = "dashed", size = 1.7) +
  geom_line(aes(y = Diffusion, color = "Diffusion"), linetype = "dotted", size = 1.7) +
  labs(title = "Black-76 Implicit Volatility",
       x = "Strike", y = "IV") +
  ylim(0.05, 0.25) +
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