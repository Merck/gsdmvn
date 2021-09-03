library(gsdmvn)
context("Testing errbeta for use in gs_design_npe()")

test_that("errbeta works does not work for K = 2 example", {
  x <- 600
  a <- errbeta(x = x, K = 2, beta = 0, theta = .1, theta1 = .1, info = 1:2, info1 = 1:2, info0 = 1:2,
               binding = FALSE, Zupper = gs_b, Zlower = gs_b, upar = 3:2, lpar = -(2:1),
               test_upper = TRUE, test_lower = TRUE, r = 18, tol = 1E-6)
  xx <-  gs_power_npe(theta = .1, theta1 = .1, info = x * (1:2), info1 = x*(1:2), info0 = x*(1:2),
                      binding = FALSE, upper = gs_b, lower = gs_b, upar = 3:2, lpar = -(2:1),
                      test_upper = TRUE, test_lower = TRUE, r = 18, tol = 1E-6) %>%
         filter(Bound == "Upper") %>%
         summarize(Power = last(Probability))
  expect_equal(1 - xx$Power, a)
})

test_that("errbeta does not work with uniroot for K = 2",{
  upper <- 600
  lower <- 500
  xx <- uniroot(errbeta, lower = lower, upper = upper,
                K = 2, beta = .1, theta = .1, theta1 = .1, info = 1:2, info1 = 1:2, info0 = 1:2,
                binding = FALSE, Zupper = gs_b, Zlower = gs_b, upar = 3:2, lpar = -(2:1),
                test_upper = TRUE, test_lower = TRUE, r = 18, tol = 1E-6)
  expect_gt(xx$root, 0)
})
