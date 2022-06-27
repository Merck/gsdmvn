test_that("Testing gridptsRcpp() vs gridpts(); easy case",{
  new <- gridptsRcpp(r = 18, mu = 0, a = -Inf, b = Inf)
  ref <- gridpts(r = 18, mu = 0, a = -Inf, b = Inf)
  expect_equal(new$w,ref$w)
  expect_equal(new$z,ref$z)
})

test_that("Testing gridptsRcpp() vs gridpts(); extreme case 1",{
  new <- gridptsRcpp(r = 18, mu = 6, a = -2, b = 0)
  ref <- gridpts(r = 18, mu = 6, a = -2, b = 0)
  expect_equal(new$w,ref$w)
  expect_equal(new$z,ref$z)
})

test_that("Testing h1Rcpp() vs h1()",{
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= rep(-.5,2), b=c(1,3))
  expect_equal(
    abs(sum(h1Rcpp(r = 18, theta = 2, I = 1, a = -Inf, b= -.5)$h) - x$lower$prob[1]),
    abs(h1(r = 18, theta = 2, I = 1, a = -Inf, b= -.5) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$lower$prob[1])
  )
  expect_equal(
    abs(sum(h1Rcpp(r = 18, theta = 2, I = 1, a = 1, b = Inf)$h) - x$upper$prob[1]),
    abs(h1(r = 18, theta = 2, I = 1, a = 1, b = Inf) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$upper$prob[1])
  )
})

test_that("Testing h1Rcpp() vs h1() for extreme case produces near 0 probability",{
  expect_equal(sum(h1Rcpp(r = 18, theta = -30, I = 1, a = -2, b = Inf)$h),
               h1(r = 18, theta = -30, I = 1, a = -2, b = Inf) %>% summarize(pupper = sum(h)) %>% as.numeric())
})

test_that("Testing hupdateRcpp() vs hupdate()", {
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= c(-.5,0), b=c(3,2))
  expect_equal(abs(x$lower$prob[2] -
                  sum(hupdateRcpp(r = 18, theta = 2, I = 2, a = -Inf, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1Rcpp(r = 18, theta = 2, I = 1, a = -.5, b = 3))$h)),
            abs(x$lower$prob[2] -
                  hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3)) %>%
                  summarize(plower2 = sum(h)) %>% as.numeric())
  )
  expect_equal(abs(x$upper$prob[2] -
                     sum(hupdateRcpp(r = 18, theta = 2, I = 2, a = 2, b = Inf, thetam1 = 2, Im1 = 1, gm1 = h1Rcpp(r = 18, theta = 2, I = 1, a = -.5, b = 3))$h)),
               abs(x$upper$prob[2] -
                  # Compare second upper crossing
                  hupdate(theta = 2, I = 2, a = 2, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3)) %>%
                  summarize(plower2 = sum(h)) %>% as.numeric()),
            1e-6
  )
})

test_that("Testing hupdateRcpp() vs hupdate() for very extreme tail cases",{
  expect_equal(sum(hupdateRcpp(r = 18, theta = 2, I = 2, a = -Inf, b = 0,
                               thetam1 = 2, Im1 = 1,
                               gm1 = h1Rcpp(r = 18, theta = -30, I = 1, a = -.5, b = 3))$h),
               hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = -30, I = 1, a = -.5, b = 3)) %>%
              summarize(plower2 = sum(h)) %>% as.numeric(),
            1e-100)
  expect_equal(sum(hupdateRcpp(r = 18, theta = 30, I = 2, a = -Inf, b = 0,
                               thetam1 = 2, Im1 = 1,
                               gm1 = h1Rcpp(r = 18, theta = 2, I = 1, a = -.5, b = 3))$h),
               hupdate(theta = 30, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3)) %>%
              summarize(plower2 = sum(h)) %>% as.numeric(),
            1e-100)
})

test_that("Testing hupdateRcpp() vs hupdate() for somewhat extreme case with non-zero probability",{
  # Test case that should have essentially 0 chance of stopping at IA1
  # and mean -sqrt(2) * 3, variance 1 at analysis 2; computing upper tail for z=1
  expect_equal(sum(hupdateRcpp(r = 18, theta = -3, I = 2, a = 0, b = Inf,
                               thetam1 = -3, Im1 = 1,
                               gm1 = h1Rcpp(r = 18, theta = -3, I = 1, a = -Inf, b = Inf))$h),
               hupdate(theta = -3, I = 2, a = 0, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf)) %>%
                 summarize(plower2 = sum(h)) %>% as.numeric(),
               tolerance=1e-100)
  # Slightly less extreme upper tail for z = -1
  # result = .000592; this is where we want to be accurate!
  expect_equal(sum(hupdateRcpp(r = 18, theta = -3, I = 2, a = -1, b = Inf,
                               thetam1 = -3, Im1 = 1,
                               gm1 = h1Rcpp(r = 18, theta = -3, I = 1, a = -Inf, b = Inf))$h),
               hupdate(theta = -3, I = 2, a = -1, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf)) %>%
                 summarize(plower2 = sum(h)) %>% as.numeric(),
               tolerance=1e-100)
})