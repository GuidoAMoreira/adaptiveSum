#' Test function
#'
#' Selects the faster method for sum approximation.
#' @export
smart_sum <- function(lFun, params = numeric(), eps = 1e-15, maxIter = 1e5,
                      logL = -Inf, n0 = 0){

  if (logL < log(0.5)) naive_sum(lFun, params, eps, maxIter, n0)
  else adapt_sum(lFun, params, eps, maxIter, logL, n0)

}
