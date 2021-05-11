#' Test function
#'
#' Wraps the adapt_sum C function for testing
#' @export
adapt_sum <- function(lFun, params = numeric(), eps = 1e-15, maxIter = 1e5,
                      logL = -Inf, n0 = 0){

  maxIter = as.integer(maxIter)
  n0 = as.integer(n0)

  if (!length(lFun)) return(list(sum = 0, n = 0))

  if (is.character(lFun)){
    if (lFun == "negbin_marginal")
      out = .Call("adapt_sum_callPrecomp",
                  1L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "noObs")
        out = .Call("adapt_sum_callPrecomp",
                    2L, params, eps, maxIter, logL, n0,
                    PACKAGE = "adaptiveSum")
    else if (lFun == "COMP")
      out = .Call("adapt_sum_callPrecomp",
                  3L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "dR0")
      out = .Call("adapt_sum_callPrecomp",
                  4L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "PL_diff")
      out = .Call("adapt_sum_callPrecomp",
                  5L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "negbin_sentinel")
      out = .Call("adapt_sum_callPrecomp",
                  6L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "poisson_sentinel")
      out = .Call("adapt_sum_callPrecomp",
                  7L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "weird_series_constL")
      out = .Call("adapt_sum_callPrecomp",
                  8L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "weird_series")
      out = .Call("adapt_sum_callPrecomp",
                  9L, params, eps, maxIter, logL, n0,
                  PACKAGE = "adaptiveSum")
  } else if(is.function(lFun)) {
    f <- function(k, Theta) lFun(k, Theta)

    out = .Call("adapt_sum",
                body(f), params, eps,
                maxIter, logL,
                n0, new.env(),
                PACKAGE = "adaptiveSum")
  } else {
    warning("Argument lFun must either be the name of a precompiled function or a function.")
    return(list(sum = 0, n = 0))
  }

  list(sum = out[[1]],
       n = out[[2]])

}
