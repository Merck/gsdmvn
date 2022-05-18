library(gsdmvn)
library(gsDesign)
library(dplyr)
context("Developer testing of gs_design_npe() inputs and fixed design")

# Parameters used repeatedly
K <- 3
info <- c(.45, .8, 1)
info0 <- info
info1 <- info
theta <- c(.2, .3, .4)
alpha <- .02
beta <- .15
upper <- gs_b
lower <- gs_b

test_that("Fixed design call to gs_design_npe() fails",{
  # Lachin book p 71 difference of proportions example
  pc <- .28 # Control response rate
  pe <- .40 # Experimental response rate
  p0 <- (pc + pe) / 2 # Ave response rate under H0
  # Information per increment of 1 in sample size
  info0 <- 1 / (p0 * (1 - p0) * 4)
  info1 <- 1 / (pc * (1 - pc) * 2 + pe * (1 - pe) * 2)
  # Result should round up to next even number = 652
  # Divide information needed under H1 by information per patient added
  n <- gs_design_npe(theta = pe - pc, info = info1, info0 = info0)$info[1] / info1
  n <- ceiling(n/2) * 2
  expect_equal(n, 652)
})

test_that("Text input for alpha in gs_design_npe() does not produce error",{
  expect_error(gs_design_npe(alpha = "this is text; should produce error"))
})

test_that("length(alpha) > 1 in gs_design_npe() does not produce error",{
  expect_error(gs_design_npe(alpha = (1:2)/20))
})

test_that("Text input for beta in gs_design_npe() does not produce error",{
  expect_error(gs_design_npe(beta = "this is text; should produce error"))
})

test_that("length(beta) > 1 does not produce error",{
  expect_error(gs_design_npe(beta = (1:2)/10))
})

test_that("alpha or beta out of range does not produce gs_design_npe() error",{
  expect_error(gs_design_npe(alpha = 2))
  expect_error(gs_design_npe(beta = -2))
  expect_error(gs_design_npe(alpha = .5, beta = .55 ))
})

test_that("info, info0, info1 must be positive increasing vectors",{
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = "text is bad"))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = -(1:3)))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = 3:1))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info0 = "text is bad"))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info0 = -(1:3)))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info0 = 3:1))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info1 = "text is bad"))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info1 = -(1:3)))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, info1 = 3:1))
})

test_that("Incorrect length for info0 or info1 in gs_design_npe() does not produce error",{
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info0 = 1:2))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info1 = 1:2))
})

test_that("theta, theta0, theta1 must be numeric vectors",{
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = "text is bad", info = info))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = 1:2, info = info))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, theta1 = c("text","is", "bad"), info = info))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, theta1 = 1:2, info = info))
})

test_that("Final effect size <= 0 in call to gs_design_npe() does not produce error",{
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = -(1:3), info = info))
})

test_that("test_upper or test_lower mis-specified does not produce error", {
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, test_upper = FALSE))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, test_upper = 1))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, test_upper = rep(TRUE, 2)))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, test_lower = 1))
  expect_error(gs_design_npe(alpha = alpha, beta = beta, theta = theta, info = info, test_lower = rep(TRUE, 2)))
})
