context("Test the adaptive sum routine")

###### Conway-Maxwell Poisson normalising constant 
COMP_lpdf <- function(k, theta){
  lambda <- theta[1]
  nu <- theta[2]
  return(
    k * log(lambda) - nu*lfactorial(k)  
  )
}

test_that("Conway-Maxwell with lambda = 20 and nu = 1", {
  Lambda <- 20
  Nu <- 1
  Theta <- c(Lambda, Nu)
  Eps <- 1E-16
  M <- 2E5
  TrueValue <- Lambda  
  logAns <- adaptiveSum::adapt_sum(lFun = COMP_lpdf, params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Conway-Maxwell with lambda = 20 and nu = 2", {
  Lambda <- 20
  Nu <- 2
  Theta <- c(Lambda, Nu)
  Eps <- 1E-16
  M <- 2E5
  TrueValue <- log(besselI(2*sqrt(Lambda), nu = 0))
  logAns <- adaptiveSum::adapt_sum(lFun = COMP_lpdf, params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Conway-Maxwell with lambda = 2 and nu = 1, higher tolerance", {
  Lambda <- 2
  Nu <- 1
  Theta <- c(Lambda, Nu)
  Eps <- 1E-10
  M <- 2E5
  TrueValue <- Lambda  
  logAns <- adaptiveSum::adapt_sum(lFun = COMP_lpdf, params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Conway-Maxwell with lambda = 2 and nu = 2, higher tolerance", {
  Lambda <- 2
  Nu <- 2
  Theta <- c(Lambda, Nu)
  Eps <- 1E-10
  M <- 2E5
  TrueValue <- log(besselI(2*sqrt(Lambda), nu = 0))
  logAns <- adaptiveSum::adapt_sum(lFun = COMP_lpdf, params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

###### Negative binomial + binomial error model  

negativeBinomial_marginalised_lpmf <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  x <- theta[4]
  ans <- ifelse(k < x,  -Inf, 
                dnbinom(k, mu = mu, size = phi, log = TRUE) + dbinom(x = x, size = k, prob = eta, log = TRUE))
  return(ans)
}

test_that("Marginalised negative binomial", {
  Mu <- 50
  Phi <- .5
  Eta <- .08
  obsX <- 20
  Theta <- c(Mu, Phi, Eta, obsX)
  Eps <- 1E-16
  TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_marginalised_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Marginalised negative binomial, larger overdispersion", {
  Mu <- 50
  Phi <- .05
  Eta <- .08
  obsX <- 20
  Theta <- c(Mu, Phi, Eta, obsX)
  Eps <- 1E-16
  TrueValue <- dnbinom(x = obsX, mu = Eta*Mu, size = Phi, log = TRUE)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_marginalised_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

###### Probability of not observing clusters, sentinel model.

probability_noObs_lpmf <- function(k, theta){
  eta <- theta[1]
  return(
    k * log(1-eta)  
  )
}

test_that("Probability of observing a zero, sentinel model lower tolerance", {
  Eta <- .08
  Theta <- Eta
  Eps <- 1E-16
  TrueValue <- -log(Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = probability_noObs_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Probability of observing a zero, sentinel model slightly higher tolerance", {
  Eta <- .08
  Theta <- Eta
  Eps <- 1E-15
  TrueValue <- -log(Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = probability_noObs_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Probability of observing a zero, sentinel model higher tolerance", {
  Eta <- .08
  Theta <- Eta
  Eps <- 1E-10
  TrueValue <- -log(Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = probability_noObs_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = -Inf, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

###### Stuttering chain marginalisation 
get_marginalised_prob <- function(r, w, p, log = FALSE, verbose = FALSE){
  ## Get Pr(S = 0 | r, w, p)
  a <- r/w
  b <- r/(w * (1 + r/w)^(w + 1)) * (1-p)
  solve_for_u <- function(b, w){
    ## We know b < 1 because a and omega are positive
    ## it suffices to find the positive root smaller than 1/omega, which
    ## is the positive root of the objective function obj_fun
    b_of_u <- function(u) u/(1 + u)^{w + 1}
    obj_fun <- function(cand) (b-b_of_u(cand))^2
    u <- optimise(obj_fun, c(0, 1/w), tol = 1e-25)$minimum
    return(u)
  }
  ustar <- solve_for_u(b, w = w)
  if(verbose) cat("u* is:", ustar, "\n")
  ans <- log1p(a)-log(a) + log(ustar) - log1p(ustar)
  if(!log) ans <- exp(ans)
  return(ans)
}
dR0 <- function(y, r, w, log = FALSE){
  if(y == 0){
    dens <- -Inf
  }else{
    c1 <- lgamma(w*y + y-1)
    c2 <- lgamma(w*y)
    c3 <- lgamma( y+1)
    c4 <- (y-1) * (log(r) - log(w))
    c5 <- (w*y + y -1) * log(1 + r/w)
    dens <- c1 - (c2+c3) + (c4-c5)
  }
  if(!log){
    dens <- exp(dens)
  }
  return(dens)
}
dR0 <- Vectorize(dR0)

R0_cluster_marginal_detection_lpmf <- function(k, theta){
  R0 <- theta[1]
  omega <- theta[2]
  x <- theta[3]
  p <- theta[4]
  ans <- ifelse(k < x, -Inf,
                dR0(y = k, r = R0, w = omega, log = TRUE) + dbinom(x = x, size = k, prob = p, log = TRUE))
  return(ans)
}

test_that("Stuttering chain marginal probability, low tolerance", {
  R0 <- .82
  Omega <- .05
  detProb <- .1
  obsX <- 0
  Theta <- c(R0, Omega, obsX, detProb)
  Eps <- 1E-20
  lgL <- log(R0) + log1p(-detProb) + (1 + Omega) * ( log1p(Omega) - log(R0 + Omega) )
  TrueValue <- get_marginalised_prob(r = R0, w = Omega, p = detProb, log = TRUE) 
  logAns <- adaptiveSum::adapt_sum(lFun = R0_cluster_marginal_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Stuttering chain marginal probability, lower tolerance", {
  R0 <- .82
  Omega <- .05
  detProb <- .1
  obsX <- 0
  Theta <- c(R0, Omega, obsX, detProb)
  Eps <- 1E-10
  lgL <- log(R0) + log1p(-detProb) + (1 + Omega) * ( log1p(Omega) - log(R0 + Omega) )
  TrueValue <- get_marginalised_prob(r = R0, w = Omega, p = detProb, log = TRUE) 
  logAns <- adaptiveSum::adapt_sum(lFun = R0_cluster_marginal_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

###### Negative binomial + sentinel model

negativeBinomial_sentinel_detection_lpmf <- function(k, theta){
  mu <- theta[1]
  phi <- theta[2]
  eta <- theta[3]
  return(
    dnbinom(k, mu = mu, size = phi, log = TRUE) + k * log1p(-eta)
  )
}

test_that("Negative binomial + sentinel observation error, lower tolerance", {
  Mu <- 20
  Phi <- 1.5
  Eta <- .08
  Theta <- c(Mu, Phi, Eta)
  Eps <- 1E-16
  TrueValue <- Phi * (log(Phi) - log(Eta*Mu + Phi)) 
  lgL <- log(Mu)- log_sum_exp(c(log(Mu), log(Phi))) + log1p(-Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_sentinel_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Negative binomial + sentinel observation error, higher tolerance", {
  Mu <- 20
  Phi <- 1.5
  Eta <- .08
  Theta <- c(Mu, Phi, Eta)
  Eps <- 1E-10
  TrueValue <- Phi * (log(Phi) - log(Eta*Mu + Phi)) 
  lgL <- log(Mu)- log_sum_exp(c(log(Mu), log(Phi))) + log1p(-Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_sentinel_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})


test_that("Negative binomial + sentinel observation error, large overdispersion, lower tolerance", {
  Mu <- 20
  Phi <- 0.5
  Eta <- .08
  Theta <- c(Mu, Phi, Eta)
  Eps <- 1E-16
  TrueValue <- Phi * (log(Phi) - log(Eta*Mu + Phi)) 
  lgL <- log(Mu)- log_sum_exp(c(log(Mu), log(Phi))) + log1p(-Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_sentinel_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

test_that("Negative binomial + sentinel observation error, large overdispersion, higher tolerance", {
  Mu <- 20
  Phi <- .5
  Eta <- .08
  Theta <- c(Mu, Phi, Eta)
  Eps <- 1E-10
  TrueValue <- Phi * (log(Phi) - log(Eta*Mu + Phi)) 
  lgL <- log(Mu)- log_sum_exp(c(log(Mu), log(Phi))) + log1p(-Eta)
  logAns <- adaptiveSum::adapt_sum(lFun = negativeBinomial_sentinel_detection_lpmf,
                                   params = Theta, eps = Eps, maxIter = 1E5, logL = lgL, n0 = 0)
  error <- robust_difference(logAns$sum, TrueValue)
  expect_true(abs(error) <= Eps)
})

