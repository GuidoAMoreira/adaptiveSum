#' Log of sum of exponentials
#'
#' @param x a vector with logarithms of the values whose sum will be computed.
#'
#' @return the log of the sum of the exponentials of the elements in x.
#' @importFrom  matrixStats logSumExp
#' @export log_sum_exp
#'
log_sum_exp <- function(x){
  ans <- matrixStats::logSumExp(lx = x)
  return(ans)
}
#' Log of 1-x
#'
#' @param x a vector of values between 0 and 1 for which to compute log(1 - x).
#'
#' @return log(1 - x).
#' @export log1m
#'
#' @examples
log1m <- function(x) {
  return(log1p(-x))
}
#' Log of 1-exp(x)
#'
#' @param a a vector of values a such that exp(a) < 1
#'
#' @return log(1-exp(a))
#' @export 
#'
log1m_exp <- function(a) {
  if (a > -0.693147) {
    return(log(-expm1(a))) 
  } else {
    return(log1m(exp(a)));
  }
}
#' Log of exp(x)-exp(y)
#'
#' @param x a scalar
#' @param y a scalar
#'
#' @return log(exp(x)-exp(y))
#' @export log_diff_exp
#'
#' @examples
log_diff_exp <- function (x, y){
  return(x + log1m_exp(y - x))
}  
#' Returns the difference between exp(x) and exp(y)
#'
#' @param x a scalar.
#' @param y a scalar.
#' @param log logical. If \code{TRUE}, the logarithmic of the difference is returned. Default is \code{FALSE}.
#'
#' @return log( max(exp(x), exp(y)) - min(exp(x), exp(y)) )
#' @export robust_difference
#'
robust_difference <- function(x, y, log = FALSE){
  sgn <- sign(x-y)
  if(log){
    ans <- log_diff_exp(max(x, y), min(x, y))
  }else{
    ans <- sgn * exp(log_diff_exp(max(x, y), min(x, y)))
  }
  return(ans)
}
#' Relative difference
#'
#' @param x a scalar.
#' @param y a scalar.
#'
#' @return relative difference, as measured by \code{all.equal} with \code{tolerance = 0}.
#' @export relative_difference
#'
relative_difference <- function(x, y){
  obj <- all.equal(as.numeric(x), as.numeric(y), tolerance = 0)
  ans <- suppressWarnings(as.numeric(gsub("Mean relative difference: ", "", obj)))
  if(is.na(ans)) ans <- 0
  return(ans)
}