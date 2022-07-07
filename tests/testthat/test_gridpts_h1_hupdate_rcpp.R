tolerance = 1e-20

test_that("Testing gridpts() vs gridpts_(); easy case",{
  new <- gridpts(r = 18, mu = 0, a = -Inf, b = Inf)
  ref <- gridpts_(r = 18, mu = 0, a = -Inf, b = Inf)
  expect_equal(new$w,ref$w, tolerance = tolerance)
  expect_equal(new$z,ref$z, tolerance = tolerance)
})

test_that("Testing gridpts() vs gridpts_(); extreme case 1",{
  new <- gridpts(r = 18, mu = 6, a = -2, b = 0)
  ref <- gridpts_(r = 18, mu = 6, a = -2, b = 0)
  expect_equal(new$w, ref$w, tolerance = tolerance)
  expect_equal(new$z, ref$z, tolerance = tolerance)
})

test_that("Testing h1() vs h1_()",{
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= rep(-.5,2), b=c(1,3))
  expect_equal(
    abs(sum(h1(r = 18, theta = 2, I = 1, a = -Inf, b= -.5)$h) - x$lower$prob[1]),
    abs(h1_(r = 18, theta = 2, I = 1, a = -Inf, b= -.5) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$lower$prob[1]), 
    tolerance = tolerance
  )
  expect_equal(
    abs(sum(h1(r = 18, theta = 2, I = 1, a = 1, b = Inf)$h) - x$upper$prob[1]),
    abs(h1_(r = 18, theta = 2, I = 1, a = 1, b = Inf) %>% summarize(plower1 = sum(h)) %>% as.numeric() - x$upper$prob[1]),
    tolerance = tolerance
  )
})

test_that("Testing h1() vs h1_() for extreme case produces near 0 probability",{
  expect_equal(sum(h1(r = 18, theta = -30, I = 1, a = -2, b = Inf)$h),
               h1_(r = 18, theta = -30, I = 1, a = -2, b = Inf) %>% summarize(pupper = sum(h)) %>% as.numeric(),
               tolerance = tolerance)
})

test_that("Testing hupdate() vs hupdate_()", {
  x <- gsDesign::gsProbability(k=2, n.I=1:2, theta = 2, a= c(-.5,0), b=c(3,2))
  expect_equal(abs(x$lower$prob[2] -
                  sum(hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h)),
               abs(x$lower$prob[2] -
                  hupdate_(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1, gm1 = h1_(theta = 2, I = 1, a = -.5, b = 3)) %>%
                  summarize(plower2 = sum(h)) %>% as.numeric()),
               tolerance = tolerance
  )
  expect_equal(abs(x$upper$prob[2] -
                     sum(hupdate(theta = 2, I = 2, a = 2, thetam1 = 2, Im1 = 1, gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h)),
               abs(x$upper$prob[2] -
                  # Compare second upper crossing
                  hupdate_(theta = 2, I = 2, a = 2, thetam1 = 2, Im1 = 1, gm1 = h1_(theta = 2, I = 1, a = -.5, b = 3)) %>%
                  summarize(plower2 = sum(h)) %>% as.numeric()),
               tolerance = tolerance
  )
})

test_that("Testing hupdate() vs hupdate_() for very extreme tail cases",{
  expect_equal(sum(hupdate(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1,
                               gm1 = h1(theta = -30, I = 1, a = -.5, b = 3))$h),
               hupdate_(theta = 2, I = 2, b = 0, thetam1 = 2, Im1 = 1,
                       gm1 = h1_(theta = -30, I = 1, a = -.5, b = 3)) %>%
              summarize(plower2 = sum(h)) %>% as.numeric(),
              tolerance = tolerance)
  expect_equal(sum(hupdate(theta = 30, I = 2, b = 0, thetam1 = 2, Im1 = 1,
                               gm1 = h1(theta = 2, I = 1, a = -.5, b = 3))$h),
               hupdate_(theta = 30, I = 2, b = 0, thetam1 = 2, Im1 = 1,
                       gm1 = h1_(theta = 2, I = 1, a = -.5, b = 3)) %>%
              summarize(plower2 = sum(h)) %>% as.numeric(),
              tolerance = tolerance)
})

test_that("Testing hupdate() vs hupdate_() for somewhat extreme case with non-zero probability",{
  # Test case that should have essentially 0 chance of stopping at IA1
  # and mean -sqrt(2) * 3, variance 1 at analysis 2; computing upper tail for z=1
  expect_equal(sum(hupdate(theta = -3, I = 2, a = 0, b = Inf,
                               thetam1 = -3, Im1 = 1,
                               gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf))$h),
               hupdate_(theta = -3, I = 2, a = 0, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1_(theta = -3, I = 1, a = -Inf, b = Inf)) %>%
                 summarize(plower2 = sum(h)) %>% as.numeric(),
               tolerance = tolerance)
  # Slightly less extreme upper tail for z = -1
  # result = .000592; this is where we want to be accurate!
  expect_equal(sum(hupdate(r = 18, theta = -3, I = 2, a = -1, b = Inf,
                               thetam1 = -3, Im1 = 1,
                               gm1 = h1(theta = -3, I = 1, a = -Inf, b = Inf))$h),
               hupdate_(theta = -3, I = 2, a = -1, b = Inf,
                       thetam1 = -3, Im1 = 1,
                       gm1 = h1_(theta = -3, I = 1, a = -Inf, b = Inf)) %>%
                 summarize(plower2 = sum(h)) %>% as.numeric(),
               tolerance = tolerance)
})
