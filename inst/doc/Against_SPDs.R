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

## ----find_spd, out.width="100%"-----------------------------------------------
spd <- FindSummedProbabilityDistribution(
  calendar_age_range_BP = c(1000, 4500), 
  rc_determinations = armit$c14_age, 
  rc_sigmas = armit$c14_sig, 
  F14C_inputs = FALSE, 
  calibration_curve = intcal20,
  plot_output = TRUE)

## ----two_normals_spd, echo = 2,  out.width="100%"-----------------------------
oldpar <- par(no.readonly = TRUE)

two_normals_spd <- FindSummedProbabilityDistribution(
  calendar_age_range_BP=c(2500, 7000),
  rc_determinations= two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  plot_output = TRUE)

par(new = TRUE, 
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(two_normals_spd$calendar_age_BP))
ylim <- c(0, 3 * max(two_normals_spd$probability))
plot(
    NULL,
    NULL,
    type = "n",
    ylim = ylim,
    xlim = xlim,
    axes = FALSE,
    xlab = NA,
    ylab = NA,
    xaxs = "i",
    yaxs = "i")
# Show true underlying calendar age density
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2

# Find and plot true exact density
truedens <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedens(
  x,
  w = weights_true,
  truemean = cluster_means_true_calBP,
  trueprec = cluster_precisions_true),
      from = 2500, to = 7000, n = 401,
      lwd = 2,
      col = "red", add = TRUE)

# Reset plotting parameters
par(oldpar)

## ----two_normals_compare, echo = FALSE,  out.width="100%"---------------------
oldpar <- par(no.readonly = TRUE)

polya_urn_output <- PolyaUrnBivarDirichlet(
  rc_determinations = two_normals$c14_age,
  rc_sigmas = two_normals$c14_sig,
  calibration_curve=intcal20,
  n_iter = 1e4,
  show_progress = FALSE)

two_normals_DPMM <- PlotPredictiveCalendarAgeDensity(
  output_data = polya_urn_output,
  show_SPD = TRUE)

par(new = TRUE, 
    mgp = c(3, 0.7, 0),
    xaxs = "i",
    yaxs = "i",
    mar = c(5, 4.5, 4, 2) + 0.1,
    las = 1)
xlim <- rev(range(two_normals_DPMM[[1]]$calendar_age_BP))
ylim <- c(0, 3 * max(two_normals_DPMM[[1]]$density_mean))
plot(
    NULL,
    NULL,
    type = "n",
    ylim = ylim,
    xlim = xlim,
    axes = FALSE,
    xlab = NA,
    ylab = NA,
    xaxs = "i",
    yaxs = "i")
# Show true underlying calendar age density
weights_true <- c(0.45, 0.55)
cluster_means_true_calBP <- c(3500, 5000)
cluster_precisions_true <- 1 / c(200, 100)^2

# Find and plot true exact density
truedens <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedens(
  x,
  w = weights_true,
  truemean = cluster_means_true_calBP,
  trueprec = cluster_precisions_true),
      from = 2500, to = 7000, n = 401,
      lwd = 2,
      col = "red", add = TRUE)

# Reset plotting parameters
par(oldpar)

