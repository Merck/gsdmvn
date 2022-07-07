
testthat::test_that("compare gridpts_ results with gsDesign::normalGrid results", {
  test1 <- gsdmvn:::gridpts_(r = 18, mu = 4, a = -Inf, b = Inf)
  x <- gsDesign::normalGrid(r = 18, bounds = c(-40, 40), mu = 4, sigma = 1)
  expect_equal(test1$w, x$gridwgts)
  expect_equal(test1$z, x$z)

  test2 <- gsdmvn::gridpts_(r = 18, mu = 2, a = -Inf, b = Inf)
  y <- gsDesign::normalGrid(r = 18, bounds = c(-40, 40), mu = 2, sigma = 1)
  expect_equal(test2$w,y$gridwgts)
  expect_equal(test2$z,y$z)
})


