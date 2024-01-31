## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dev='png',
  fig.width=7, 
  fig.height=5.5 # ,
  # dev.args=list(antialias = "none")
)

## ----setup--------------------------------------------------------------------
library(carbondate)

## ----calculate_gr_multiple, results=FALSE, eval=FALSE-------------------------
#  all_outputs <- list()
#  for (i in 1:3) {
#    set.seed(i + 1)
#    all_outputs[[i]] <- PolyaUrnBivarDirichlet(
#        kerr$c14_age, kerr$c14_sig, intcal20, n_iter = 1e4)
#  }
#  PlotGelmanRubinDiagnosticMultiChain(all_outputs)

## ----out.width= "100%", echo = FALSE------------------------------------------
knitr::include_graphics("figures-convergence/calculate_gr_multiple-1.png")

## ----calculate_gr, results=FALSE, eval=FALSE----------------------------------
#  set.seed(3)
#  output <- PolyaUrnBivarDirichlet(
#      kerr$c14_age, kerr$c14_sig, intcal20, n_iter = 2e4)
#  
#  PlotGelmanRubinDiagnosticSingleChain(output, n_burn = 5e3)

## ----out.width= "100%", echo = FALSE------------------------------------------
knitr::include_graphics("figures-convergence/calculate_gr-1.png")

## ----calculate_polya_kerr, results=FALSE, eval=FALSE--------------------------
#  outputs <- list()
#  for (i in 1:3) {
#    set.seed(i+1)
#    outputs[[i]] <- PolyaUrnBivarDirichlet(
#        rc_determinations = kerr$c14_age,
#        rc_sigmas = kerr$c14_sig,
#        calibration_curve=intcal20,
#        n_iter = 1e4)
#    outputs[[i]]$label <- paste("Seed =", i)
#  }
#  PlotPredictiveCalendarAgeDensity(
#    outputs, n_posterior_samples = 500, denscale = 2, interval_width = "1sigma")

## ----out.width= "100%", echo = FALSE------------------------------------------
knitr::include_graphics("figures-convergence/calculate_polya_kerr-1.png")

## ----calculate_polya_normals, results=FALSE, eval=FALSE-----------------------
#  outputs <- list()
#  for (i in 1:3) {
#    set.seed(i + 1)
#    outputs[[i]] <- PolyaUrnBivarDirichlet(
#    rc_determinations = two_normals$c14_age,
#    rc_sigmas = two_normals$c14_sig,
#    calibration_curve=intcal20,
#    n_iter = 1e4)
#    outputs[[i]]$label <- paste("Seed =", i)
#  }
#  PlotPredictiveCalendarAgeDensity(
#    outputs, n_posterior_samples = 500, denscale = 2, interval_width = "1sigma")

## ----out.width= "100%", echo = FALSE------------------------------------------
knitr::include_graphics("figures-convergence/calculate_polya_normals-1.png")

## ----calculate_kld, results=FALSE, eval=FALSE---------------------------------
#  set.seed(50)
#  output <- WalkerBivarDirichlet(
#      rc_determinations = kerr$c14_age,
#      rc_sigmas = kerr$c14_sig,
#      calibration_curve=intcal20,
#      n_iter = 1e5)
#  
#  PlotConvergenceData(output)

## ----out.width= "100%", echo = FALSE------------------------------------------
knitr::include_graphics("figures-convergence/calculate_kld-1.png")

