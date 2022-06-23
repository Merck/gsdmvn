test_that("gridpts fails comparison with gsDesign package; easy case",{
  new <- gridptsRcpp(r = 18, mu = 0, a = -Inf, b = Inf)
  ref <- gridpts(r = 18, mu = 0, a = -Inf, b = Inf)
  expect_equal(new$w,ref$w)
  expect_equal(new$z,ref$z)
})

test_that("gridpts fails comparison with gsDesign package; extreme case 1",{
  new <- gridptsRcpp(r = 18, mu = 6, a = -2, b = 0)
  ref <- gridpts(r = 18, mu = 6, a = -2, b = 0)
  expect_equal(new$w,ref$w)
  expect_equal(new$z,ref$z)
})
