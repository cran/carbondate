#' Model Occurrence of Multiple Radiocarbon Samples as a
#' Variable-Rate Poisson Process
#'
#' @description
#' This function calibrates a set of radiocarbon (\eqn{{}^{14}}C) samples, and provides an estimate
#' of how the underlying rate at which the samples occurred varies over calendar time (including any
#' specific changepoints in the rate). The method can be used as an alternative approach to summarise
#' calendar age information contained in a set of related \eqn{{}^{14}}C samples, enabling inference
#' on the latent \emph{activity} rate which led to their creation.
#'
#' It takes as an input a set of \eqn{{}^{14}}C determinations and associated \eqn{1\sigma}
#' uncertainties, as well as the radiocarbon age calibration curve to be used. The (calendar age)
#' occurrence of these radiocarbon (\eqn{{}^{14}}C) samples is modelled as a Poisson process. The
#' underlying rate of this Poisson process \eqn{\lambda(t)}, which represents the rate at which the
#' samples occurred, is considered unknown and assumed to vary over calendar time.
#'
#' Specifically, the sample occurrence rate \eqn{\lambda(t)} is modelled as piecewise constant, but with
#' an unknown number of changepoints, which can occur at unknown times. The value of \eqn{\lambda(t)}
#' between any set of changepoints is also unknown. The function jointly calibrates the given \eqn{{}^{14}}C
#' samples under this model, and simultaneously provides an estimate of \eqn{\lambda(t)}. Fitting is performed
#' using reversible-jump MCMC within Gibbs.
#'
#' It returns estimates for the calendar age of each individual radiocarbon sample in the input set;
#' and broader output on the estimated value of \eqn{\lambda(t)} which can be used by other library functions.
#' Analysis of the estimated changepoints in the piecewise \eqn{\lambda(t)} permits wider inference on whether
#' the occurrence rate of samples changed significantly at any particular calendar time and, if so, when and
#' by how much.
#'
#' For more information read the vignette: \cr
#' \code{vignette("Poisson-process-modelling", package = "carbondate")}
#'
#' @inheritParams WalkerBivarDirichlet
#'
#' @param calendar_age_range (Optional) Overall minimum and maximum calendar ages (in cal yr BP) permitted
#' for the set of radiocarbon samples, i.e., `calendar_age_range[1]` < \eqn{\theta_1, \ldots, \theta_n} <
#' `calendar_age_range[1]`. This is used to bound the start and end of the Poisson process (so no events
#' will be permitted to occur outside this interval). If not selected then automated selection will be
#' made based on given `rc_determinations` and value of `bounding_range_prob_cutoff`
#'
#' @param calendar_grid_resolution The spacing of the calendar age grid on which to restrict
#' the potential calendar ages of the samples, e.g., calendar ages of samples are limited to being one of
#' `t, t + resolution, t + 2 * resolution, ...` Default is 1 (i.e., all calendar years in the overall
#' calendar range are considered). Primarily used to speed-up code if have large range, when may select
#' larger resolution.
#'
#' @param prior_h_shape,prior_h_rate (Optional) Prior for the value of the Poisson Process rate (the height `rate_h`)
#' in any specific interval:
#' \deqn{\textrm{rate}\_\textrm{h} \sim \textrm{Gamma}(
#' \textrm{shape} = \textrm{prior}\_\textrm{h}\_\textrm{shape},
#' \textrm{rate} = \textrm{prior}\_\textrm{h}\_\textrm{rate}).}
#' If they are both `NA` then `prior_h_shape` is selected to be 1 (so `rate_h` follows an Exponential
#' distribution) and `prior_h_rate` chosen adaptively (internally) to match `n_observations`.
#'
#' @param prior_n_internal_changepoints_lambda Prior mean for the number of internal changepoints
#' in the rate \eqn{\lambda(t)}.
#' \deqn{\textrm{n}\_\textrm{internal}\_\textrm{changepoints} \sim
#' \textrm{Po}(\textrm{prior}\_\textrm{n}\_\textrm{internal}\_\textrm{changepoints}\_\textrm{lambda})}
#'
#' @param k_max_internal_changepoints Maximum permitted number of internal changepoints
#'
#' @param rescale_factor_rev_jump Factor weighting probability of dimension change
#' in the reversible jump update step for Poisson process `rate_h` and `rate_s`
#'
#' @param bounding_range_prob_cutoff Probability cut-off for choosing the bounds for the
#' potential calendar ages for the observations
#'
#' @param initial_n_internal_changepoints Number of internal changepoints to initialise MCMC sampler with.
#' The default is 10 (so initialise with diffuse state). Will place these initial changepoints uniformly at random
#' within overall calendar age range.
#'
#' @param grid_extension_factor If you adaptively select the `calendar_age_range` from `rc_determinations`, how
#' far you wish to extend the grid beyond this adaptive minimum and maximum. The final range will be extended
#' (equally on both sides) to cover `(1 + grid_extension_factor) * (calendar_age_range)`
#'
#' @param use_fast,fast_approx_prob_cutoff A flag to allow trimming the calendar age likelihood (i.e., reducing the
#' range of calendar ages) for each individual sample to speed up the sampler. If `TRUE` (default), for each
#' individual sample, those tail calendar ages (in the overall `calendar_age_grid`) with very small likelihoods will
#' be trimmed (speeding up the updating of the calendar ages). If `TRUE` the probability cut-off to remove the tails
#' is `fast_approx_prob_cutoff`.
#'
#' @return A list with 7 items. The first 4 items contain output of the model, each of
#' which has one dimension of size \eqn{n_{\textrm{out}} =
#' \textrm{floor}( n_{\textrm{iter}}/n_{\textrm{thin}}) + 1}. The rows in these items store
#' the state of the MCMC from every \eqn{n_{\textrm{thin}}}\eqn{{}^\textrm{th}} iteration:
#'
#' \describe{
#'  \item{`rate_s`}{A list of length \eqn{n_{\textrm{out}}} each entry giving the current set of
#'  (calendar age) changepoint locations in the piecewise-constant rate \eqn{\lambda(t)}.}
#'  \item{`rate_h`}{A list of length \eqn{n_{\textrm{out}}} each entry giving the current set of
#'  heights (values for the rate) in each piecewise-constant section of \eqn{\lambda(t)}.}
#'  \item{`calendar_ages`}{An \eqn{n_{\textrm{out}}} by \eqn{n_{\textrm{obs}}}
#'     matrix. Gives the current estimate for the calendar age of each individual
#'     observation.}
#'  \item{`n_internal_changes`}{A vector of length \eqn{n_{\textrm{out}}} giving the current number of
#'  internal changes in the value of \eqn{\lambda(t)}.}
#' }
#' where \eqn{n_{\textrm{obs}}} is the number of radiocarbon observations, i.e.,
#' the length of `rc_determinations`.
#'
#' The remaining items give information about input data, input parameters (or
#' those calculated) and update_type
#'
#' \describe{
#'  \item{`update_type`}{A string that always has the value "RJPP".}
#'  \item{`input_data`}{A list containing the \eqn{{}^{14}}C data used, and the name of
#'  the calibration curve used.}
#'  \item{`input_parameters`}{A list containing the values of the fixed
#'  parameters `pp_cal_age_range`, `prior_n_internal_changepoints_lambda`,
#'  `k_max_internal_changepoints`, `prior_h_shape`, `prior_h_rate`, `rescale_factor_rev_jump`,
#'   `calendar_age_grid`, `calendar_grid_resolution`, `n_iter` and `n_thin`.}
#' }
#' @export
#'
#' @examples
#' # NOTE: This example is shown with a small n_iter to speed up execution.
#' # When you run ensure n_iter gives convergence (try function default).
#'
#' pp_output <- PPcalibrate(
#'     pp_uniform_phase$c14_age,
#'     pp_uniform_phase$c14_sig,
#'     intcal20,
#'     n_iter = 100,
#'     show_progress = FALSE)
PPcalibrate <- function(
    rc_determinations,
    rc_sigmas,
    calibration_curve,
    F14C_inputs = FALSE,
    n_iter = 1e5,
    n_thin = 10,
    use_F14C_space = TRUE,
    show_progress = TRUE,
    calendar_age_range = NA,
    calendar_grid_resolution = 1,
    prior_h_shape = NA,
    prior_h_rate = NA,
    prior_n_internal_changepoints_lambda = 3,
    k_max_internal_changepoints = 30,
    rescale_factor_rev_jump = 0.9,
    bounding_range_prob_cutoff = 0.001,
    initial_n_internal_changepoints = 10,
    grid_extension_factor = 0.1,
    use_fast = TRUE,
    fast_approx_prob_cutoff = 0.001) {

  arg_check <- .InitializeErrorList()
  .CheckInputData(arg_check, rc_determinations, rc_sigmas, F14C_inputs)
  .CheckCalibrationCurve(arg_check, calibration_curve, NA)
  .CheckIterationParameters(arg_check, n_iter, n_thin)
  .CheckFlag(arg_check, F14C_inputs)
  .CheckFlag(arg_check, use_F14C_space)
  .CheckFlag(arg_check, show_progress)
  if (!any(is.na(calendar_age_range))) { .CheckNumberVector(arg_check, calendar_age_range, len = 2) }
  .CheckNumber(arg_check, calendar_grid_resolution, lower = 0)
  .CheckPriorHShapeAndPriorHRate(arg_check, prior_h_shape, prior_h_rate)
  .CheckInteger(arg_check, k_max_internal_changepoints, lower = 1)
  .CheckNumber(arg_check, bounding_range_prob_cutoff, lower = 0, upper = 0.01)
  .CheckInteger(
    arg_check, initial_n_internal_changepoints, upper = k_max_internal_changepoints - 1)
  .CheckNumber(arg_check, grid_extension_factor, lower = 0)
  .CheckFlag(arg_check, use_fast)
  .CheckNumber(arg_check, fast_approx_prob_cutoff, lower = 0, upper = 0.01)
  .ReportErrors(arg_check)

  ##############################################################################
  # Save input data
  input_data <- list(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    F14C_inputs = F14C_inputs,
    calibration_curve_name = deparse(substitute(calibration_curve)))

  ##############################################################################
  # Convert the scale of the initial determinations to F14C or C14_age as appropriate
  # if they aren't already
  if (F14C_inputs == use_F14C_space) {
    rc_determinations <- as.double(rc_determinations)
    rc_sigmas <- as.double(rc_sigmas)
  } else if (F14C_inputs == FALSE) {
    converted <- .Convert14CageToF14c(rc_determinations, rc_sigmas)
    rc_determinations <- converted$f14c
    rc_sigmas <- converted$f14c_sig
  } else {
    converted <- .ConvertF14cTo14Cage(rc_determinations, rc_sigmas)
    rc_determinations <- converted$c14_age
    rc_sigmas <- converted$c14_sig
  }

  # Find initial calendar_age_range for Poisson Process
  # Have converted rc_determinations to use_F14C_space already
  bounds_calendar_range <- .FindBoundingCalendarRange(
    rc_determinations = rc_determinations,
    rc_sigmas = rc_sigmas,
    calibration_curve = calibration_curve,
    F14C_inputs = use_F14C_space,
    prob_cutoff = bounding_range_prob_cutoff)

  if (any(is.na(calendar_age_range))) {
    min_potential_calendar_age <- min(bounds_calendar_range)
    max_potential_calendar_age <- max(bounds_calendar_range)

    # Extend by grid_extension_factor
    mean_potential_calendar_age <- mean(bounds_calendar_range)
    min_potential_calendar_age <- (
      min_potential_calendar_age - grid_extension_factor * (
        mean_potential_calendar_age - min_potential_calendar_age
      )
    )
    max_potential_calendar_age <- (
      max_potential_calendar_age + grid_extension_factor * (
        max_potential_calendar_age - mean_potential_calendar_age
      )
    )

    # Extend to discretised calendar_grid_resolution
    min_potential_calendar_age <- (
      floor(min_potential_calendar_age/calendar_grid_resolution) *
        calendar_grid_resolution
    )
    max_potential_calendar_age <- (
      ceiling(max_potential_calendar_age/calendar_grid_resolution) *
        calendar_grid_resolution
    )

    # Do not extend beyond calibration curve calendar limits
    min_potential_calendar_age <- max(min_potential_calendar_age,
                                      min(calibration_curve$calendar_age_BP))
    max_potential_calendar_age <- min(max_potential_calendar_age,
                                      max(calibration_curve$calendar_age_BP))
  } else {
    ## Create calendar_age_grid covering potential calendar ages
    min_potential_calendar_age <- min(calendar_age_range)
    max_potential_calendar_age <- max(calendar_age_range)
    # Provide warning if these are narrower than the estimated bounds
    if(min_potential_calendar_age > min(bounds_calendar_range)) {
      warning("The minimum of your selected calendar age range may not cover the age range of the samples. Consider reducing the minimum age range." , immediate. = TRUE)
    }
    if(max_potential_calendar_age < max(bounds_calendar_range)) {
      warning("The maximum of your selected calendar age range may not cover the age range of the samples. Consider increasing the maximum age range." , immediate. = TRUE)
    }
  }

  calendar_age_grid <- seq(
    min_potential_calendar_age,
    max_potential_calendar_age,
    by = calendar_grid_resolution) # May end before max_potential_calendar_age

  # Ensure end of calendar_age_grid extends at least to max_potential_calendar_age
  # If not extend calendar_age_grid and adjust max_potential_calendar_age to match
  if((max(calendar_age_grid) != max_potential_calendar_age)) {
    if( !any(is.na(calendar_age_range)) ) { # Only inform if have user selected values
      warning("Extending calendar age range for Poisson Process as selected grid resolution doesn't divide age range.", immediate. = TRUE)
    }
    max_potential_calendar_age <- (
      calendar_age_grid[length(calendar_age_grid)] + calendar_grid_resolution
    )
    cc <- c(calendar_age_grid, max_potential_calendar_age)
  }

  calendar_age_interval_length <- max_potential_calendar_age - min_potential_calendar_age
  n_determinations <- length(rc_determinations)

  initial_estimate_mean_rate <- n_determinations / calendar_age_interval_length

  if(is.na(prior_h_shape)) {
    # Choose exponential distribution
    # and match mean with initial_estimate_mean_rate above
    prior_h_shape <- 1 # This is equivalent to an exponential distribution
    prior_h_rate <- 1 / initial_estimate_mean_rate
  }

  #############
  ## Create initial change points and heights for Poisson process rate
  ## Randomly place initial_n_internal_changepoints
  initial_rate_s <- sort(
    c(
      min_potential_calendar_age,
      stats::runif(initial_n_internal_changepoints,
            min = min_potential_calendar_age,
            max = max_potential_calendar_age),
      max_potential_calendar_age
    )
  )

  # Create diverse initial_rate_h centred on initial_estimate_mean_rate
  # Requires safe setting as as will crash if any sampled rate_h values are zero
  initial_rate_h <- rep(initial_estimate_mean_rate,
                        times = initial_n_internal_changepoints + 1)
  initial_rate_h_log_multiplier <- stats::runif(n = initial_n_internal_changepoints + 1,
                                                min = -1,
                                                max = 1)
  initial_rate_h <- exp(initial_rate_h_log_multiplier) * initial_rate_h

  initial_integrated_rate <- .FindIntegral(
    rate_s = initial_rate_s,
    rate_h = initial_rate_h
  )


  ##############################################################################
  # Save input parameters
  input_parameters <- list(
    pp_cal_age_range = c(min_potential_calendar_age,
                         max_potential_calendar_age),
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    prior_h_shape = prior_h_shape,
    prior_h_rate = prior_h_rate,
    rescale_factor_rev_jump = rescale_factor_rev_jump,
    calendar_age_grid = calendar_age_grid,
    calendar_grid_resolution = calendar_grid_resolution,
    n_iter = n_iter,
    n_thin = n_thin)

  ####################################
  ## Create matrix of calendar_likelihoods (stored in main as not updated throughout samples)
  ## Different from .ProbabilitiesForSingleDetermination as not normalised
  ## and can pass theta on different grid to calibration_curve
  likelihood_calendar_ages_from_calibration_curve <- mapply(
    .CalendarAgeLikelihoodGivenCurve,
    rc_determinations,
    rc_sigmas,
    MoreArgs = list(
      theta = calendar_age_grid,
      F14C_inputs = use_F14C_space,
      calibration_curve = calibration_curve))

  if (use_fast) {
    likelihood_values <- list()
    likelihood_offsets <- c()
    for (i in 1:dim(likelihood_calendar_ages_from_calibration_curve)[2]) {
      ret <- .FindTrimmedVectorAndIndices(
        likelihood_calendar_ages_from_calibration_curve[, i], fast_approx_prob_cutoff)
      likelihood_values[[i]] <- ret$values
      likelihood_offsets[i] <- as.integer(ret$offset - 1)
    }
  } else  {
    likelihood_values <- NA
    likelihood_offsets <- NA
  }

  num_observations <- length(rc_determinations)

  # Set starting values to be initialised ones
  rate_s <- initial_rate_s
  rate_h <- initial_rate_h
  integrated_rate <- initial_integrated_rate

  prob_move <- .FindMoveProbability(
    prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
    k_max_internal_changepoints = k_max_internal_changepoints,
    rescale_factor = rescale_factor_rev_jump)

  ####################################
  # Create storage for output
  n_out <- floor(n_iter / n_thin) + 1

  rate_s_out <- list(rate_s)
  rate_h_out <- list(rate_h)
  n_internal_changes <- rep(NA, length = n_out)
  theta_out <- matrix(NA, nrow = n_out, ncol = num_observations)

  output_index <- 1
  n_internal_changes[output_index] <- length(rate_h) - 1

  ## Store calendar_ages given initial_rate_s and initial_rate_h
  ## (sample from exactly using Gibbs)
  calendar_ages <- .UpdateCalendarAges(
    likelihood_calendar_ages_from_calibration_curve,
    likelihood_values,
    likelihood_offsets,
    calendar_age_grid,
    rate_s,
    rate_h,
    use_fast)
  theta_out[output_index, ] <- calendar_ages


  #####################################
  # Perform MCMC - RJMCMC within Gibbs
  # Consist of iterating between:
  #    i) Updating calendar_ages given lambda (rate_s, rate_h) and rc_determinations
  #    ii) Updating lambda (rate_s, rate_h) given calendar_ages using RJ MCMC

  if (show_progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)
  }


  for(iter in 1:n_iter) {

    ## Step 1: Update calendar_ages given rate_s and rate_h (sample from exactly using Gibbs)
    calendar_ages <- .UpdateCalendarAges(
      likelihood_calendar_ages_from_calibration_curve,
      likelihood_values,
      likelihood_offsets,
      calendar_age_grid,
      rate_s,
      rate_h,
      use_fast)

    ## Step 2: Update rate_s and rate_h given current calendar_ages (using RJMCMC)
    updated_poisson_process <- .UpdatePoissonProcessRateRevJump(
      theta = calendar_ages,
      rate_s = rate_s,
      rate_h = rate_h,
      integrated_rate = integrated_rate,
      prior_h_shape = prior_h_shape,
      prior_h_rate = prior_h_rate,
      prior_n_internal_changepoints_lambda = prior_n_internal_changepoints_lambda,
      prob_move = prob_move
    )

    rate_s <- updated_poisson_process$rate_s
    rate_h <- updated_poisson_process$rate_h
    integrated_rate <- updated_poisson_process$integrated_rate

    # Store output
    if (iter %% n_thin == 0) {
      output_index <- output_index + 1
      n_internal_changes[output_index] <- length(rate_h) - 1
      rate_h_out[[output_index]] <- rate_h
      rate_s_out[[output_index]] <- rate_s
      theta_out[output_index, ] <- calendar_ages
    }

    if (show_progress) {
      if (iter %% 100 == 0) {
        utils::setTxtProgressBar(progress_bar, iter)
      }
    }
  }

  return_list <- list(
    rate_s = rate_s_out,
    rate_h = rate_h_out,
    calendar_ages = theta_out,
    n_internal_changes = n_internal_changes,
    update_type = "RJPP",
    input_data = input_data,
    input_parameters = input_parameters)

  if (show_progress) close(progress_bar)
  return(return_list)
}


.UpdateCalendarAges <- function(
  likelihood_calendar_ages_from_calibration_curve,
  likelihood_values,
  likelihood_offsets,
  calendar_age_grid,
  rate_s,
  rate_h,
  use_fast) {
  if (use_fast) {
    calendar_ages <- .TrimmedUpdateCalendarAgesGibbs(
      calendar_age_grid = calendar_age_grid,
      rate_s = rate_s,
      rate_h = rate_h,
      likelihood_values = likelihood_values,
      likelihood_offsets = likelihood_offsets)
  } else {
    calendar_ages <- .UpdateCalendarAgesGibbs(
      likelihood_calendar_ages_from_calibration_curve = likelihood_calendar_ages_from_calibration_curve,
      calendar_age_grid = calendar_age_grid,
      rate_s = rate_s,
      rate_h = rate_h
    )
  }
  return(calendar_ages)
}
