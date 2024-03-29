test_that("WalkerBivarDirichlet gives expected result - C14 space and same parameters", {
  # Here we load the predictive density calculated from the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare the predictive density using the Kullback-Leibler
  # divergence. From experimenting around we expect a divergence of less the 1e-3 to say
  # that the results match.
  skip_on_cran()
  load(test_path("fixtures", "walker_input.rda"))
  set.seed(seednum)
  walker_output <- WalkerBivarDirichlet(
    rc_determinations=x,
    rc_sigmas=xsig,
    F14C_inputs = FALSE,
    calibration_curve=intcal20,
    use_F14C_space = FALSE,
    lambda=lambda,
    nu1=nu1,
    nu2=nu2,
    A=A,
    B=B,
    alpha_shape=cprshape,
    alpha_rate=cprrate,
    n_iter=niter,
    n_thin=nthin,
    calendar_ages=inittheta,
    slice_width=slicew,
    slice_multiplier=m,
    n_clust=kstar,
    sensible_initialisation = FALSE,
    mu_phi = stats::median(inittheta))

  load(test_path("fixtures", "walker_output.rda"))

  # Here tempx is the calendar age sequence used to find the predictive density in the legacy code.
  pred_dens <- FindPredictiveCalendarAgeDensity(walker_output, tempx, 5000)

  expect_lt(.KLD(postden, pred_dens$density_mean), 1e-3)

})

test_that("WalkerBivarDirichlet gives expected result - sensible init, F14C space", {
  # Here we load input and output data from the functions in the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare results to test that the function is
  # working correctly
  skip_on_cran()
  load(test_path("fixtures", "walker_input.rda"))

  set.seed(seednum)

  walker_output <- WalkerBivarDirichlet(
    rc_determinations=x,
    rc_sigmas=xsig,
    F14C_inputs = FALSE,
    use_F14C_space = TRUE,
    calibration_curve=intcal20,
    n_iter=1e5,
    n_thin=10)

  load(test_path("fixtures", "walker_output.rda"))

  # Here tempx is the calendar age sequence used to find the predictive density in the legacy code.
  pred_dens <- FindPredictiveCalendarAgeDensity(walker_output, tempx, 5000)

  expect_lt(.KLD(postden, pred_dens$density_mean), 1e-3)
})

test_that("WalkerBivarDirichlet gives expected result - sensible init, F14C space, F14C inputs", {
  # Here we load input and output data from the functions in the original
  # nonparametric calibration code (that we know is working correctly). We run
  # the function again here and compare results to test that the function is
  # working correctly
  skip_on_cran()
  load(test_path("fixtures", "walker_input.rda"))

  set.seed(seednum)

  f14c_determinations <- exp(-x / 8033)
  f14c_sig <- xsig * f14c_determinations / 8033

  walker_output <- WalkerBivarDirichlet(
    rc_determinations=f14c_determinations,
    rc_sigmas=f14c_sig,
    F14C_inputs = TRUE,
    use_F14C_space = TRUE,
    calibration_curve=intcal20,
    n_iter=1e5,
    n_thin=10)

  load(test_path("fixtures", "walker_output.rda"))

  # Here tempx is the calendar age sequence used to find the predictive density in the legacy code.
  pred_dens <- FindPredictiveCalendarAgeDensity(walker_output, tempx, 5000)

  expect_lt(.KLD(postden, pred_dens$density_mean), 1e-3)
})


test_that("WalkerBivarDirichlet gives error if no sensible init and not all params", {
  # This should fail with multiple assertions. Just test one of them.
  expect_error(
    WalkerBivarDirichlet(
      rc_determinations=kerr$c14_age,
      rc_sigmas=kerr$c14_sig,
      F14C_inputs = FALSE,
      calibration_curve=intcal20,
      sensible_initialisation = FALSE),
    regexp = "lambda must be a number"
  )
})

