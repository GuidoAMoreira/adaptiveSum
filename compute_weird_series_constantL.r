source("R/truncation_package.r")
source("R/truncation_aux.r")
library(adaptiveSum)

weird_series <- function(k, theta){
  s <- theta[1]
  if(k==0) return(-Inf)
  return(
    -(2*log(k) + k*log(s) )
  )
}
weird_series <- Vectorize(weird_series)

S <- 1.5
Theta <- S
Eps <- 1E-14
TrueValue <- log(copula::polylog(z = 1/S, s = 2))
lgL <- -log(S)

result <- compare_approximations(weird_series, theta = Theta,
                                 exact = TrueValue, epsilon = Eps, max_iter = 1E5, logL = lgL,
                                 compiled_function_name = "weird_series_constL")

result

library(microbenchmark)
mit <- 1E5

timing <- microbenchmark(
  naive = approx_naive(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  naive_threshold = approx_naive_under(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  doubling = approx_doubling(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  adaptive = approx_adaptive(weird_series, theta = Theta, epsilon = Eps, max_iter = mit),
  adaptiveC = adapt_sum(weird_series, params = Theta, eps = Eps, maxIter = mit),
  adaptiveComp = adapt_sum("weird_series_constL", params = Theta, eps = Eps, maxIter = mit)
)

result
timing
