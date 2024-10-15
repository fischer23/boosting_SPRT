# boosting_SPRT

This repository contains the R code for reproducing the results in the paper "Improving the (approximate) sequential probability ratio test by avoiding overshoot". In the aforementioned paper, we introduce a general approach to improve test martingale based sequential tests. In particular, we uniformly improve the power-one SPRT with threshold 1/\alpha and propose an alternative to Wald's approximated two-sided thresholds that often needs less samples while providing provable type I error control. 

Files:

boosting_functions.R:          Contains the functions for boosting martingale factors. 

calc_boosting_factors.R:       Calculates the boosting factors reported in the Sections 3.1, 4.3 and 5 of the paper. 

simple.R:                      Generates the simulated data for Section 3.1 and saves it in the file "simple.rda" in the results folder.

comp_alt.R:                    Generates the simulated data for Section 3.2 and saves it in the file "comp_alt.rda" in the results folder.

functions_CS.R:                Contains auxiliary functions for the computation of boosted confidence sequences.

CS.R:                          Generates the simulated data for Section 4.1 and saves it in the file "CS.rda" in the results folder.

WoR.R:                         Generates the simulated data for Section 4.2 and saves it in the file "WoR.rda" in the results folder.

futility.R:                    Generates the simulated data for Section 5 and saves it in the file "futility.rda" in the results folder.

plot_generator.R:               Uses the simulated data in the results folder to create all plots of the paper and save them in the results folder.
