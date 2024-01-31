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

## ----illustrate_set, out.width="100%", echo = FALSE, fig.cap = paste("_Illustration of our challenge. We observe a set of related samples, each with a ^14^C determination (shown here as the red ticks on the radiocarbon age axis). Can we jointly calibrate the samples, and summarise the combined calendar age information that they provide?_")----
set.seed(13)
oldpar <- par(no.readonly = TRUE)

par(
  mgp = c(3, 0.7, 0),
  xaxs = "i",
  yaxs = "i",
  mar = c(5, 4.5, 4, 2) + 0.1,
  las = 1)

calcurve <- intcal20

# Sample 14C determinations between artificial start and end dates 
mincal <- 2000
maxcal <- 3000
nsamp <- 50
theta_true <- sample(mincal:maxcal, size = nsamp, replace = TRUE)

# Sample some 14C determinations
xsig <- rep(25, nsamp)
x <- rnorm(nsamp, mean = approx(calcurve$calendar_age_BP, calcurve$c14_age, theta_true)$y,
           sd = sqrt(xsig^2 + (approx(calcurve$calendar_age_BP, calcurve$c14_sig, theta_true)$y)^2) )

plot(calcurve$calendar_age_BP, calcurve$c14_age, col = "blue", 
     ylim = range(x) + c(-4,4) * max(xsig), xlim = c(maxcal, mincal) + c(100, -100),
     xlab = "Calendar Age (cal yr BP)", 
     ylab = expression(paste("Radiocarbon age (", ""^14, "C ", "yr BP)")),
     type = "l", main = expression(paste(""^14,"C Calibration and Summarisation")))
calcurve$ub <- calcurve$c14_age + 2 * calcurve$c14_sig
calcurve$lb <- calcurve$c14_age - 2 * calcurve$c14_sig
lines(calcurve$calendar_age_BP, calcurve$ub, lty = 2, col = "blue" )
lines(calcurve$calendar_age_BP, calcurve$lb, lty = 2, col = "blue")
polygon(c(rev(calcurve$calendar_age_BP), calcurve$calendar_age_BP), c(rev(calcurve$lb), calcurve$ub), col=rgb(0,0,1,.3), border=NA)
rug(x, side = 2, ticksize = 0.03, lwd = 1, col = "red")
legend_labels <- c(
    substitute(paste(""^14, "C determination ")),
    "IntCal20",
    expression(paste("2", sigma, " interval")))
lty <- c(1, 1, 2)
lwd <- c(1, 1, 1)
pch <- c(NA, NA, NA)
col <- c(grDevices::rgb(1, 0, 0, .5), "blue", "blue")
legend(
    "topright", legend = legend_labels, lty = lty, lwd=lwd, pch = pch, col = col)

# Reset plotting parameters
par(oldpar)


