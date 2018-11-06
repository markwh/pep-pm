# stats.R
# Evaluation statistics

RRMSE <- function (pred, meas) 
  sqrt(mean((pred - meas)^2/meas^2))

NRMSE <- function (pred, meas) 
  sqrt(mean((meas - pred)^2))/mean(meas)


geomMean <- function(x) {
  exp(mean(log(x)))
}
