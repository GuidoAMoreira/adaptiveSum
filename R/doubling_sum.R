#' Test function
#'
#' Wraps the naive_sum C function for testing
#' @export
doubling_sum <- function(lFun, params = numeric(), eps = 1e-15, N_start = 20,
                         c = 2, maxIter = 1e5, n0 = 0){

  maxIter = as.integer(maxIter)
  n0 = as.integer(n0)
  N_start = as.integer(N_start)
  c = as.integer(c)

  if (!length(lFun)) return(list(sum = 0, n = 0))

  if (is.character(lFun)){
    if (lFun == "negbin_marginal")
      out = .Call("doubling_sum_callPrecomp",
                  1L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "noObs")
      out = .Call("doubling_sum_callPrecomp",
                  2L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "COMP")
      out = .Call("doubling_sum_callPrecomp",
                  3L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "dR0")
      out = .Call("doubling_sum_callPrecomp",
                  4L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "PL_diff")
      out = .Call("doubling_sum_callPrecomp",
                  5L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "negbin_sentinel")
      out = .Call("doubling_sum_callPrecomp",
                  6L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "poisson_sentinel")
      out = .Call("doubling_sum_callPrecomp",
                  7L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "weird_series_constL")
      out = .Call("doubling_sum_callPrecomp",
                  8L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
    else if (lFun == "weird_series")
      out = .Call("doubling_sum_callPrecomp",
                  9L, params, eps, N_start, c, maxIter, n0,
                  PACKAGE = "adaptiveSum")
  } else if(is.function(lFun)) {
    f <- function(k, Theta) lFun(k, Theta)

    out = .Call("doubling_sum",
                body(f), params, eps, N_start, c,
                maxIter,
                n0, new.env(),
                PACKAGE = "adaptiveSum")
  } else {
    warning("Argument lFun must either be the name of a precompiled function or a function.")
    return(list(sum = 0, n = 0))
  }

  list(sum = out[[1]],
       n = out[[2]])

}
