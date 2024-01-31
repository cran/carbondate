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

## ----illustrate_PP, fig.cap = paste("_Modelling sample occurence using a Poisson process with a latent occurence/activity rate.  Left Panel: An illustration of a Poisson process with two changepoints in the occurrence rate. The samples (shown as a rug in purple) occur at random calendar times, but proportional to the underlying occurrence rate. There are therefore more samples between (2700, 2300) cal yr BP  than other times. Right Panel: Each sample has a ^14^C age (shown as black ticks on the y-axis). We wish to reconstruct the underlying Poisson process rate given only these ^14^C determination._"), out.width="100%", fig.width=10, fig.height=4, echo = FALSE----

oldpar <- par(no.readonly = TRUE)

par(
  mfrow = c(1,2),
  mgp = c(3, 0.7, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(5, 4.5, 4, 2) + 0.1,
  las = 1)

set.seed(16)

# Artificial start and end dates 
mincal <- 1950
maxcal <- 3100
t <- seq(mincal, maxcal, by = 1)

# Create a simulated lambda (piecewise constant)
lambda <- c(rep(0.1, 350), rep(0.45, 400), rep(0.2, length(t) - (400 + 350)))

# Choose beta(mutiplier)
beta <- 0.2
truerate <- beta * lambda

# Sample calendar ages (according to Poisson process)
nsamp <- rpois(1, sum(truerate))
truetheta <- sample(t, nsamp, replace = TRUE, prob = truerate)

# Create plot showing rate and calendar ages
plot(c(min(t), t, max(t)), 
     c(0,lambda,0), 
     lty = 1, lwd = 2, 
     col = "red", type = "l", 
     ylim = c(0, max(lambda) + 0.1), 
     xlim = c(maxcal, mincal) + c(100, -100), yaxs = "i",
     xlab = "Calendar Age (cal yr BP)", 
     ylab = expression(
       paste("Occurrence Rate ", lambda, "(", theta, ")" )),
     main = "Occurence of Samples")
rug(truetheta, side = 1, 
    col = "purple", 
    lwd = 2, ticksize = 0.05,
    quiet = TRUE)
abline(v = t[350], lty = 2, lwd = 2)
abline(v = t[350 + 400], lty = 2, lwd = 2)


# Sample some 14C determinations
rc_sigmas <- rep(25, nsamp)
rc_determinations <- rnorm(nsamp, 
                           mean = approx(intcal20$calendar_age_BP, 
                                         intcal20$c14_age, 
                                         truetheta)$y,
                           sd = sqrt(rc_sigmas^2 + (approx(intcal20$calendar_age_BP, 
                                                           intcal20$c14_sig, 
                                                           truetheta)$y)^2) 
           )

# Plot 14C determinations against IntCal20
plot(intcal20$calendar_age_BP, 
     intcal20$c14_age, 
     col = "blue", 
     ylim = range(rc_determinations) + c(-4,4) * max(rc_sigmas), 
     xlim = c(maxcal, mincal) + c(100, -100),
     xlab = "Calendar Age (cal yr BP)", 
     ylab = expression(paste("Radiocarbon age (", ""^14, "C ", "yr BP)")),
     type = "l", main = expression(paste(""^14,"C Changepoint Analysis")))
lines(intcal20$calendar_age_BP, 
      intcal20$c14_age + 2 * intcal20$c14_sig,
      lty = 2, 
      col = "blue" )
lines(intcal20$calendar_age_BP, 
      intcal20$c14_age - 2 * intcal20$c14_sig,
      lty = 2, 
      col = "blue" )
polygon(c(rev(intcal20$calendar_age_BP), intcal20$calendar_age_BP), 
        c(rev(intcal20$c14_age - 2 * intcal20$c14_sig), 
          intcal20$c14_age + 2 * intcal20$c14_sig), 
        col=rgb(0,0,1,.3), border=NA)
rug(rc_determinations, side = 2, quiet = TRUE)

# Plot the true rate along the bottom
par(new = TRUE)
plot(c(min(t), t, max(t)), c(0,truerate,0), lty = 1, col = "red", type = "l", 
     ylim = c(0, 10 * max(truerate)), xlim = c(maxcal, mincal) + c(100, -100),
     axes = FALSE, xlab = NA, ylab = NA, yaxs = "i")

# Reset plotting parameters
par(oldpar)

## ----example, out.width="100%"------------------------------------------------
set.seed(15)

# Set initial values 
n_observed <- 40
rc_sigmas <- rep(15, n_observed)

# Create artificial rc_determinations
calendar_age_range <- c(300, 700) 
observed_age_range <- c(500, 550) 
true_theta <- seq(
      from = observed_age_range[1],
      to = observed_age_range[2],
      length = n_observed)
intcal_mean <- approx(
      x = intcal20$calendar_age_BP,
      y = intcal20$c14_age,
      xout = true_theta)$y
intcal_sd <- approx(
      x = intcal20$calendar_age_BP,
      y = intcal20$c14_sig,
      xout = true_theta)$y
rc_determinations <- rnorm(
      n = n_observed,
      mean = intcal_mean,
      sd = sqrt(rc_sigmas^2 + intcal_sd^2))

# Fit the model
PP_fit_output <- PPcalibrate(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = intcal20,
    calendar_age_range = calendar_age_range,
    calendar_grid_resolution = 1,
    n_iter = 1e5,
    n_thin = 10,
    show_progress = FALSE)

## ----meanrate, out.width="100%"-----------------------------------------------
PlotPosteriorMeanRate(PP_fit_output)

## ----changepoint_number, out.width="100%"-------------------------------------
PlotNumberOfInternalChanges(PP_fit_output)

## ----changepoint_locations, out.width="100%"----------------------------------
PlotPosteriorChangePoints(PP_fit_output)

## ----changepoint_rates, out.width="100%"--------------------------------------
PlotPosteriorHeights(PP_fit_output)

## ----PP_plot_individual, out.width="100%"-------------------------------------
PlotCalendarAgeDensityIndividualSample(
  9, PP_fit_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)

## ----changepoint_locations_AD, out.width="100%"-------------------------------
PlotPosteriorChangePoints(PP_fit_output, 
                          plot_cal_age_scale = "AD")

