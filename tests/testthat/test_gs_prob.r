library(gsdmvn)
context("Multivariate probabilities for group sequential design")

test_that("Testing gs_prob() vsresults using mvtnorm package",{

# Version of gs_prob using mvtnorm package ##########################################
gsProbx <- function(theta, upper, lower, info,
                   algorithm = mvtnorm::GenzBretz(maxpts = 50000, abseps = 1e-05)){
  K <- length(info)
  # If only a single theta argument, use gsProbability
  if (length(theta) ==1){
    d <- gsDesign::gsProbability(k = K, n.I = info, theta = theta,
                                 a = pmax(lower, -20), b = pmin(upper, 20))
    return(tibble::tibble(Analysis = rep(1:K, 2),
                          Bound = c(rep("Upper", K), rep("Lower", K)),
                          Z = c(upper, lower),
                          Probability = c(cumsum(d$upper$prob),
                                          cumsum(d$lower$prob)),
                          theta = theta) %>% filter(abs(Z) < Inf)
    )
  }
  # If multiple theta values, use pmvnorm (for now!)
  mu <- theta * sqrt(info)
  ti <- matrix(rep(info, K), nrow = K, byrow = TRUE)
  corr <- sqrt(pmin(ti, t(ti)) /  pmax(ti, t(ti)))
  upperProb <- rep(NA, K)
  lowerProb <- rep(NA, K)
  for(k in seq_along(info)){
    if(k==1){
      upperProb[1] <- if(upper[1] < Inf) {pnorm(upper[1], mean = mu[1], lower.tail = FALSE)}else{0}
      lowerProb[1] <- if(lower[1] > -Inf){pnorm(lower[1], mean = mu[1])}else{0}
    }else{
      upperProb[k] <- if(upper[k]< Inf){
        mvtnorm::pmvnorm(lower = c(lower[1:(k-1)], upper[k]),
                         upper = c(upper[1:(k-1)], Inf),
                         corr = corr[1:k,1:k],
                         mean = mu[1:k],
                         algorithm = algorithm)
      }else{0}
      lowerProb[k] <- if(lower[k] > -Inf){
        mvtnorm::pmvnorm(lower = c(lower[1:(k-1)], -Inf),
                         upper = c(upper[1:(k-1)], lower[k]),
                         corr = corr[1:k,1:k],
                         mean = mu[1:k],
                         algorithm = algorithm)
      }else{0}
    }
  }
  return(tibble::tibble(Analysis = rep(1:K, 2),
                        Bound = c(rep("Upper", K), rep("Lower", K)),
                        Z = c(upper, lower),
                        Probability = c(cumsum(upperProb),
                                        cumsum(lowerProb)),
                        theta = rep(theta, 2)) %>% filter(abs(Z) < Inf)
  )
}
###################################################################################

# Set up an example
theta = 1:3
upper = c(3,2.5,2)
lower = c(-2,-1,0)
info = 3:5

## mvtnorm version
set.seed(1239)
x <- gsProbx(theta = theta, upper = upper, lower =  lower, info = info)
y <- gs_prob(theta = theta, upar = upper, lpar =  lower, info = info)
expect_gt(1E-6, max(abs(x$Probability - y$Probability)))
})
##################################################################################


test_that("Testing gs_prob() vs gsDesign::gsProbability(); single effect size",{
  x <- gsDesign::gsProbability(k = 5, theta = 1, n.I = 1:5, b = rep(2.5,5), a = rep(0,5))
  y <- gs_prob(theta = 1, info = 1:5, upper= gs_b, upar = rep(2.5, 5), lower=gs_b, lpar=rep(0,5))$Probability
  diff  <- max(abs(c(cumsum(x$upper$prob), cumsum(x$lower$prob)) - y))
  expect_lt(diff, 1E-6)
})


test_that("Testing gs_prob() vs gsDesign::gsProbability(); efficacy only",{
  x <- gsDesign::gsProbability(k = 5, theta = 1, n.I = 1:5, b = rep(2.5,5), a = rep(-20,5))
  y <- gs_prob(theta = 1, info = 1:5, upper = gs_b, upar = rep(2.5, 5), lower = gs_b, lpar = rep(-Inf,5))$Probability
  diff  <- max(abs(cumsum(x$upper$prob) - y[1:5]))
  expect_lt(diff, 1E-6)
})

test_that("Testing gs_prob() vs gsDesign::gsProbability(); futility only",{
  x <- gsDesign::gsProbability(k = 5, theta = 1, n.I = 1:5, b = rep(20,5), a = -2:2)
  y <- gs_prob(theta = 1, info = 1:5, upper = gs_b, upar = rep(Inf, 5), lower =  gs_b, lpar = -2:2)$Probability
  diff  <- max(abs(cumsum(x$lower$prob) - y[6:10]))
  expect_lt(diff, 1E-6)
})

test_that("Testing gs_prob() vs gsDesign::gsProbability(); futility only followed by efficacy only",{
  x <- gsDesign::gsProbability(k = 5, theta = 0, n.I = 1:5, b = c(20, 20, rep(qnorm(.975),3)), a = c(rep(qnorm(.05), 2), rep(-20, 3)))
  y <- gs_prob(theta = 0, info = 1:5, upper = gs_b, upar = c(Inf, Inf, rep(qnorm(.975),3)), lower = gs_b, lpar = c(rep(qnorm(.05), 2), rep(-Inf, 3)))$Probability
  diff  <- max(abs(c(cumsum(x$upper$prob[3:5]), cumsum(x$lower$prob[1:2])) - y[3:7]))
  expect_lt(diff, 1E-6)
})

test_that("Testing gs_prob() vs pnorm() for relatively small p-value",{
  y <- gs_prob(theta = 0, info = (1:5)/5,
               upper = gs_b, upar = rep(Inf, 5),
               lower = gs_b, lpar = c(rep(-Inf, 4), -3))$Probability %>% max()
  expect_equal(y, pnorm(-3), tol = 1.07E-7)
})

test_that("Testing gs_prob() vs pnorm() for relatively large p-value",{
  y <- gs_prob(theta = 0, info = (1:5)/5,
               upper = gs_b, upar = rep(Inf, 5),
               lower = gs_b, lpar = c(rep(-Inf, 4), 3))$Probability %>% max()
  expect_equal(y, pnorm(3), tol = 1.03E-6)
})
