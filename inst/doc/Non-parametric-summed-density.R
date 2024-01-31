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
set.seed(5)

## ----illustrate_mixture, fig.cap = paste("_An illustration of building up a potentially complex distribution $f(\\theta)$ using mixtures of normals. Left Panel: A simple mixture of three (predominantly disjoint) normal clusters (blue dashed lines) results in an overall $f(\\theta)$ that is tri-modal (solid red). Right Panel: Overlapping normal clusters (blue dashed lines) can however create more complex $f(\\theta)$ distributions (solid red)._"), out.width="100%", fig.width=10, fig.height=4, echo = FALSE----

oldpar <- par(no.readonly = TRUE)

par(
  mfrow = c(1,2),
  mgp = c(3, 0.7, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(5, 4.5, 2, 2) + 0.1,
  las = 1)

# Plots summed normal distributions
t <- 1:100
wtrue <- c(0.3, 0.4, 0.3)
phitrue <- c(40, 60, 80)
tautrue <- 1/c(5, 4, 8)^2

# Find and plot true exact density
truedensv2 <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}

curve(truedensv2(x, w = wtrue, truemean = phitrue, trueprec = tautrue), 
      from = min(t), to = max(t), n = 401, 
      ylim = c(0, 0.05), 
      ylab = expression(paste("f(", theta, ")")),
      xlab = "Calendar Age (cal yr BP)",
      xlim = rev(c(min(t), max(t))),
      col = "red", lwd = 2)
# Add in constituent normals
for(i in 1:length(wtrue)) {
  curve(truedensv2(x, w = wtrue[i], truemean = phitrue[i], trueprec = tautrue[i]), 
        from = min(t), to = max(t), n = 401, add = TRUE,
        ylim = c(0, 0.06), 
        xlim = rev(c(min(t), max(t))),
        col = "blue", lty = 2)
}


## Version two - overlapping
wtrue <- c(0.2, 0.3, 0.1, 0.2, 0.2)
phitrue <- c(50, 60, 70, 40, 30)
tautrue <- 1/c(5, 4, 8, 4, 7)^2

curve(truedensv2(x, w = wtrue, truemean = phitrue, trueprec = tautrue), 
      from = min(t), to = max(t), n = 401, 
      ylim = c(0, 0.05), 
      ylab = expression(paste("f(", theta, ")")),
      xlab = "Calendar Age (cal yr BP)",
      xlim = rev(c(min(t), max(t))),
      col = "red", lwd = 2)
# Add in constituent normals
for(i in 1:length(wtrue)) {
  curve(truedensv2(x, w = wtrue[i], truemean = phitrue[i], trueprec = tautrue[i]), 
        from = min(t), to = max(t), n = 401, add = TRUE,
        ylim = c(0, 0.06), 
        xlim = rev(c(min(t), max(t))),
        col = "blue", lty = 2)
}

# Reset plotting parameters
par(oldpar)

## ----calculate_walker, results=FALSE------------------------------------------
polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = kerr$c14_age,
  rc_sigmas = kerr$c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 5)

## ----calculate_neal, results=FALSE--------------------------------------------
walker_output <- WalkerBivarDirichlet(
  rc_determinations = kerr$c14_age,
  rc_sigmas = kerr$c14_sig,
  calibration_curve = intcal20,
  n_iter = 1e5,
  n_thin = 5)

## ----plot_density, out.width="100%"-------------------------------------------
densities <- PlotPredictiveCalendarAgeDensity(
  output_data = list(polya_urn_output, walker_output),
  denscale = 2.5)

## ----plot_density_2, out.width="100%"-----------------------------------------
densities <- PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  denscale = 2.5,
  show_SPD = TRUE,
  interval_width = "bespoke",
  bespoke_probability = 0.8,
  plot_14C_age = FALSE)

## ----plot_individual, out.width="100%"----------------------------------------
PlotCalendarAgeDensityIndividualSample(21, polya_urn_output)

## ----plot_individual_with_hpd, out.width="100%"-------------------------------
PlotCalendarAgeDensityIndividualSample(
  21, polya_urn_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)

## ----plot_clusters, out.width="50%", fig.width=5, fig.height=6----------------
PlotNumberOfClusters(output_data = polya_urn_output)

## ----plot_density_2_AD, out.width="100%"--------------------------------------
densities <- PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  show_SPD = TRUE,
  plot_cal_age_scale = "AD")

