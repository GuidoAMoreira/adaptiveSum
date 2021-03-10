#' Test function
#'
#' Wraps the adapt_sum C function for testing
#' @importFrom matrixStats logSumExp
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

  n = out[[2]]
  list(sum = matrixStats::logSumExp(out[[1]][1:(n+2)]),
       n = n,
       out = out)

}
