library(testthat)

testthat::test_that("expect equal with gsDesign::gsProbability outcome for signle effect size, efficacy, futility, k=3",{


  test1 <- gs_prob(theta = 0,
                   upar = gsDesign::gsDesign(k = 3,sfu = gsDesign::sfLDOF)$upper$bound,
                   lpar = gsDesign::gsDesign(k = 3,sfl = gsDesign::sfLDOF)$lower$bound,
                   upper = gs_b,
                   lower = gs_b,
                   info = 1:3)

  x <- gsDesign::gsProbability(
    k = 3,
    a = gsDesign::gsDesign(k = 3,sfl = gsDesign::sfLDOF)$lower$bound,
    b = gsDesign::gsDesign(k = 3,sfu = gsDesign::sfLDOF)$upper$bound,
    n.I = gsDesign::gsDesign(k = 3,sfu = gsDesign::sfLDOF)$n.I,
    theta = 0
  )

  expect_equal(object = test1[1:3,]$Probability, expected = cumsum(x$upper$prob), tolerance = 0.0001)
  expect_equal(object = test1[4:6,]$Probability, expected = cumsum(x$lower$prob), tolerance = 0.0001)
})


testthat::test_that("expect equal with mvtnorm outcome at final analysis, k=2",{

  P <- function(alpha.t, alpha.ia, r, b){
    temp = mvtnorm::pmvnorm(lower = c(-Inf, b),
                            upper = c(qnorm(1 - alpha.ia), Inf),
                            corr = rbind(c(1, sqrt(r)), c(sqrt(r), 1)))
    return(alpha.t - alpha.ia - temp)
  }

  alpha.t <- 0.025
  t <- 0.5
  r <- t
  b.ia <- sfLDOF(alpha = alpha.t, t = t)
  alpha.ia <- b.ia$spend
  b <- uniroot(P, c(1.96, 4), alpha.t = alpha.t, alpha.ia = alpha.ia, r = r)
  p <- 1- pnorm(b$root)

  test2 <- gs_prob(theta = 0,
                   upar = gsDesign::gsDesign(k = 2,sfu = gsDesign::sfLDOF)$upper$bound,
                   lpar = gsDesign::gsDesign(k = 2,sfl = gsDesign::sfLDOF)$lower$bound,
                   upper = gs_b,
                   lower = gs_b,
                   info = 1:2)

  expect_equal(object = test2[2,]$Probability, expected = p, tolerance = 0.001)
})
