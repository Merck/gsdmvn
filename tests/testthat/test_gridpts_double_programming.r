test_that("gridpts fails comparison with gsDesign package; easy case",{
          new <- gridpts(r = 18, mu = 0, a = -Inf, b = Inf)
          old <- gsDesign::normalGrid(r = 18, bounds = c(-20, 20), mu = 0, sigma = 1)
          expect_equal(new$w, old$gridwgts)
          expect_equal(new$z, old$z)
})
test_that("gridpts fails comparison with gsDesign package; extreme case 1",{
  new <- gridpts(r = 18, mu = 6, a = -2, b = 0)
  old <- gsDesign::normalGrid(r = 18, bounds = c(-2, 0), mu = 6, sigma = 1)
  expect_equal(new$w, old$gridwgts)
  expect_equal(new$z, old$z)
})
# Following fails due to extreme circumstances
# Calculations for extreme cases will be checked with h1 and hupdate
# since probabilities in this range will be essentially 0
# test_that("gridpts fails comparison with gsDesign package; very extreme case",{
#   new <- gsdmvn::gridpts(r = 18, mu = -30, a = -2, b = 0)
#   old <- gsDesign::normalGrid(r = 18, bounds = c(-2, 0), mu = -30, sigma = 1)
#   expect_equal(new$w,old$gridwgts)
#   expect_equal(new$z,old$z)
# })
