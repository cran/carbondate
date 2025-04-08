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
# Run plotting function and assign output to PP_posterior_mean_rate_2sigma
PP_posterior_mean_rate_2sigma <- PlotPosteriorMeanRate(PP_fit_output)
# Note: You can change the PP line widths using the plot_lwd argument 

# Look at PP_posterior_mean_rate_2sigma (posterior mean and 2 sigma intervals) 
head(PP_posterior_mean_rate_2sigma)


## ----no_assign_meanrate, eval = FALSE-----------------------------------------
# # Run plotting function
# PlotPosteriorMeanRate(PP_fit_output)

## ----assign_meanrate, out.width="100%", fig.show='hide'-----------------------
# Run plotting function and assign output to PP_posterior_mean_rate_2sigma
PP_posterior_mean_rate_2sigma <- PlotPosteriorMeanRate(PP_fit_output)
# Will also recreate plot above (but not shown in vignette)

# Look at PP_posterior_mean_rate_2sigma (posterior mean and 2 sigma intervals) 
head(PP_posterior_mean_rate_2sigma)

## ----changepoint_number, out.width="100%"-------------------------------------
PlotNumberOfInternalChanges(PP_fit_output)

## ----changepoint_locations, out.width="100%"----------------------------------
PlotPosteriorChangePoints(PP_fit_output)
# Can add an n_changes argument, e.g., n_changes = c(2, 3, 4)
# if want to condition on a different number of changes 

## ----changepoint_rates, out.width="100%"--------------------------------------
PlotPosteriorHeights(PP_fit_output)
# As above can add an n_changes argument, e.g., n_changes = c(2, 3, 4)
# if want to condition on a different number of changes 

## ----PP_plot_individual, out.width="100%"-------------------------------------
PlotCalendarAgeDensityIndividualSample(
  9, PP_fit_output, show_hpd_ranges = TRUE, show_unmodelled_density = TRUE)

## ----changepoint_locations_AD, out.width="100%"-------------------------------
PlotPosteriorChangePoints(PP_fit_output, 
                          plot_cal_age_scale = "AD")

## ----find_meanrate, out.width="100%"------------------------------------------
# Calculating 2 sigma (95.4%) intervals on posterior mean occurrence rate 
PP_posterior_mean_rate_2sigma <- FindPosteriorMeanRate(
  PP_fit_output,
  calendar_age_sequence = seq(300, 500, by = 1),
  interval_width = "2sigma")

# Look at posterior mean with 2 sigma probability interval
head(PP_posterior_mean_rate_2sigma)

## ----plot_conditionalmeanrate, out.width="100%"-------------------------------
# Conditional on TWO internal changes in the occurrence rate, 
# Calculate and plot the posterior mean rate over time 
# (with its 2 sigma intervals)
conditional_2_changes_posterior_mean_rate <- PlotPosteriorMeanRate(
  PP_fit_output,
  n_changes = 2) # here n_changes must have length one (i.e., a single number)

# Look at conditional posterior mean with 2 sigma probability interval
head(conditional_2_changes_posterior_mean_rate)

## ----plot_ind_realisation, out.width="100%"-----------------------------------
# Choose some nice plotting colours (from Okabe-Ito)
realisation_colours <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442" )

# Plot 5 random realisations from posterior
PlotRateIndividualRealisation(
     PP_fit_output,
     n_realisations = 5,
     plot_realisations_colour = realisation_colours)

