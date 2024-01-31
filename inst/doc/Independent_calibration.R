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

## ----single_calibration, out.width="100%"-------------------------------------
calibration_result <- CalibrateSingleDetermination(
  rc_determination = 1413, 
  rc_sigma = 25, 
  F14C_inputs = FALSE, 
  calibration_curve = intcal20,
  plot_output = TRUE)

## ----single_calibration_SH, out.width="100%"----------------------------------
calibration_result <- CalibrateSingleDetermination(
  rc_determination = 1413, 
  rc_sigma = 25, 
  F14C_inputs = FALSE, 
  calibration_curve = shcal20,
  plot_output = TRUE)

## ----single_calibration_AD, out.width="100%"----------------------------------
calibration_result <- CalibrateSingleDetermination(
  rc_determination = 1413, 
  rc_sigma = 25, 
  F14C_inputs = FALSE, 
  calibration_curve = intcal20,
  plot_output = TRUE,
  plot_cal_age_scale = "AD")

