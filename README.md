# MScSOFR

This repository contains R code for the simulations and calculations behind the Master's Thesis "Term Structure Modeling of SOFR: A Theoretical Approach to Explore Distributional Variations in Scheduled Jumps" by Anna Caroline Lestrøm Bertelsen, cand.scient.oecon from University of Copenhagen.


The file LibrariesAndParameters.R contains code for the loading of packages used and parameter values.

The file OvernightRates.R contains code for Figure 1 of the overnight rates displayed in Chapter 2 with data accessed from the Federal Reserve Bank of St. Louis (FRED).

The file MixedExpApproximatesNormal.R contains code for the approximation of the normal distribution, N(\mu^P, (\sigma^P)^2), by the mixed-exponential distribution Mixed-Exp(\bm{\eta}^P, \bm{\theta}^P), as well as for the plot generated for Figure 2.

The file SkellamPMF.R contains code for Figure 3.

The file JumpModelPDF.R contains code for the the Monte-Carlo based probability density function of each of the four model specifications plotted in Figure 4.

The file Simulating1M.R contains code for the simulations plotted in Figure 5.

The file Simulating1Y.R contains code for the simulations plotted in Figure 6.

The file SkellamZCBODEs.R contains code for the ordinary differential equations used to calculate bond prices with Skellam-specified scheduled jump sizes. This file thus holds the functions for Eq. (6.16)-(6.17).

The file ZCB.R contains code for the calculation of zero-coupon bond prices for each of the four model specifications as well as for the plot displayed in Figure 7.

The file NormalXDMoments.R contains functions used to calculate the conditional moments of the Gaussian-specified scheduled jump sizes given in Eq. (5.46)-(5.48).

The file SkellamXDMoments.R contains the function used to calculate the conditional moments of the Skellam-specified scheduled jump sizes given in Eq. (7.16).

The file IVTimeNN.R contains code for the calculation of the nearest at-the-money One-Month SOFR call option price and associated Black-76 implied volatility with a model specification assuming Gaussian scheduled jump sizes and normal unscheduled jump sizes.

The file IVTimeNME.R contains code for the calculation of the nearest at-the-money One-Month SOFR call option price and associated Black-76 implied volatility with a model specification assuming Gaussian scheduled jump sizes and mixed-exponential unscheduled jump sizes.

The file IVTimeSN.R contains code for the calculation of the nearest at-the-money One-Month SOFR call option price and associated Black-76 implied volatility with a model specification assuming Skellam scheduled jump sizes and normal unscheduled jump sizes.

The file IVTimeSME.R contains code for the calculation of the nearest at-the-money One-Month SOFR call option price and associated Black-76 implied volatility with a model specification assuming Skellam scheduled jump sizes and mixed-exponential unscheduled jump sizes.

The file IVTimeDiffusion.R contains code for the calculation of the nearest at-the-money One-Month SOFR call option price and associated Black-76 implied volatility with a model specification assuming no jumps and only a pure diffusion process as specified by the spread process, $s_t$.

The file IVTimePlots.R contains code for the plots generated in Figures 8 and 9.

The file SkellamIVODEs.R contains code for the ordinary differential equations used to calculate the price of a One-Month SOFR call option expiring in one year with Skellam-specified scheduled jump sizes. This file thus holds the functions for Eq. (8.54)-(8.55).

The file IVStrikeNN.R contains code for the calculation of a One-Month SOFR call option price expiring in one year and associated Black-76 implied volatility with a model specification assuming Gaussian scheduled jump sizes and normal unscheduled jump sizes.

The file IVStrikeNME.R contains code for the calculation of a One-Month SOFR call option price expiring in one year and associated Black-76 implied volatility with a model specification assuming Gaussian scheduled jump sizes and mixed-exponential unscheduled jump sizes.

The file IVStrikeSN.R contains code for the calculation of a One-Month SOFR call option price expiring in one year and associated Black-76 implied volatility with a model specification assuming Skellam scheduled jump sizes and normal unscheduled jump sizes.

The file IVStrikeSME.R contains code for the calculation of a One-Month SOFR call option price expiring in one year and associated Black-76 implied volatility with a model specification assuming Skellam scheduled jump sizes and mixed-exponential unscheduled jump sizes.

The file IVStrikeDiffusion.R contains code for the calculation of a One-Month SOFR call option price expiring in one year and associated Black-76 implied volatility with a model specification assuming no jumps and only a pure diffusion process as specified by the spread process, $s_t$.

The file IVStrikeMC.R contains code for the Monte-Carlo calculated price of a One-Month SOFR call option expiring in one year and associated Black-76 implied volatility with a model specification assuming no jumps and only a pure diffusion process as specified by the spread process, $s_t$.

The file IVStrikePlots.R contains code for the plots generated in Figure 10.

