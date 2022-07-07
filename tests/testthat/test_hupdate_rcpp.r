library(gsdmvn)
context("Updated grid and weights for numerical integration")

test_that("Testing hupdate_() vs known results from gsProbability", {
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= c(-.5,0), b=c(3,2))
  expect_lt(abs(x$lower$prob[2] -
    sum(hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h)),
    1e-6
  )
  expect_lt(abs(x$upper$prob[2] -
                 # Compare second upper crossing
                 sum(hupdate(theta = 2, I = 2, a = 2, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h)),
            1e-6
  )
})
test_that("Testing hupdate for very extreme tail cases",{
  expect_lt(sum(hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = -30, I = 1, a = -.5, b = 3))$h),
            1e-100)
  expect_lt(sum(hupdate(theta = 30, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h),
            1e-100)
})
test_that("Testing somewhat extreme case with non-zero probability",{
  # Test case that should have essentially 0 chance of stopping at IA1
  # and mean -sqrt(2) * 3, variance 1 at analysis 2; computing upper tail for z=1
  expect_equal(sum(hupdate(theta = -3, I = 2, a = 0, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf))$h),
               pnorm(0, mean = -sqrt(2) * 3, lower.tail = FALSE),
               tolerance=1e-7)
  # Slightly less extreme upper tail for z = -1
  # result = .000592; this is where we want to be accurate!
  expect_equal(sum(hupdate(theta = -3, I = 2, a = -1, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf))$h),
               pnorm(-1, mean = -sqrt(2) * 3, lower.tail = FALSE),
               tolerance=1e-7)
})
