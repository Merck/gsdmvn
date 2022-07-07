library(gsdmvn)
context("Initial grid and weights for numerical integration")

test_that("Testing h1_() vs known results from gsProbability",{
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= rep(-.5,2), b=c(1,3))
  expect_lt(
    abs(gsdmvn:::h1_(theta = 2, I = 1, b= -.5) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$lower$prob[1]),
    1e-6
  )
  expect_lt(
    abs(gsdmvn:::h1_(theta = 2, I = 1, a = 1) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$upper$prob[1]),
    1e-6
  )
})
test_that("Testing that extreme case produces near 0 probability",{
  expect_lt(gsdmvn:::h1_(theta = -30, I = 1, a = -2, b = Inf) %>% summarize(pupper = sum(h)) %>% as.numeric(), 1e-100)
})
