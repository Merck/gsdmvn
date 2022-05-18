library(gsdmvn)
library(gsDesign)
context("Testing spending designs vs gsDesign package")

# Parameters used repeatedly
K <- 3
timing <- c(.45, .8, 1)
sfu <- sfPower
sfupar <- 4
sfl <- sfHSD
sflpar <- 2
delta <- .2
alpha <- .02
beta <- .15

test_that("One-sided design fails to reproduce gsDesign package bounds", {
  gsd <- gsDesign(
    test.type = 1, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
    delta = delta, alpha = alpha, beta = beta)
  
  gsdv <- gs_design_npe(
    theta = delta, info = timing, beta = beta,
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_b,
    lpar = rep(-Inf, K)
    ) %>% 
    filter(Bound == "Upper" & hypothesis == "H1")
  
  expect_equal(gsd$upper$bound, gsdv$Z, tolerance = 7e-6)
  expect_equal(gsd$n.I, gsdv$info, tolerance = .001)
  
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(
    theta = 0, info = gsdv$info,
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_b,
    lpar = rep(-Inf, K)
    ) %>% 
    filter(Bound == "Upper" & hypothesis == "H1")
  
  expect_equal(gsdv0$Probability, sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})


test_that("Two-sided symmetric design fails to reproduce gsDesign test.type=2 bounds", {
  
  gsd <- gsDesign(
    test.type = 2, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
    delta = delta, alpha = alpha, beta = beta, tol = 1e-6)
  
  gsdv <- gs_design_npe(
    theta = delta, info = timing, beta = beta,
    theta1 = rep(0,3), # Use this for lower bound spending under null hypothesis
    binding = FALSE, # Use this for 2-sided symmetric design
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_spending_bound,
    lpar = list(sf = sfu, total_spend = alpha, param = sfupar),
    tol = 1e-6
    ) %>% 
    filter(hypothesis == "H1")

  expect_equal(gsd$upper$bound, gsdv$Z[1:K], tolerance = 7e-6)
  expect_equal(gsd$lower$bound, gsdv$Z[(K+1):(2*K)], tolerance = 7e-6)
  expect_equal(gsd$n.I, gsdv$info[1:K], tolerance = .04) # While tolerance should not be problematic, it seems large
  
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(
    theta = 0, info = gsdv$info[1:K],
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_spending_bound,
    lpar = list(sf = sfu, total_spend = alpha, param = sfupar)
    ) %>% 
    filter(hypothesis == "H1")
  expect_equal(gsdv0$Probability[1:K], 
               sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})

test_that("Two-sided asymmetric design fails to reproduce gsDesign test.type=3 bounds", {
  
  gsd <- gsDesign(
    test.type = 3, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
    delta = delta, alpha = alpha, beta = beta)
  
  gsdv <- gs_design_npe(
    theta = delta, info = timing, beta = beta,
    binding = TRUE, # Use this for test.type=3 and 5
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_spending_bound,
    lpar = list(sf = sfl, total_spend = beta, param = sflpar)
    ) %>% 
    filter(hypothesis == "H1")

  expect_equal(gsd$upper$bound, gsdv$Z[1:K], tolerance = 7e-6)
  expect_equal(gsd$lower$bound, gsdv$Z[(K+1):(2*K)], tolerance = 9e-6)
  expect_equal(gsd$n.I, gsdv$info[1:K], tolerance = .04) # While tolerance should not be problematic, it seems large
  
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(
    theta = 0, info = gsdv$info[1:K],
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_spending_bound,
    lpar = list(sf = sfu, total_spend = alpha, param = sfupar)
    ) %>% 
    filter(hypothesis == "H1")
  
  expect_equal(gsdv0$Probability[1:K], sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})

test_that("Two-sided asymmetric design fails to reproduce gsDesign test.type=4 bounds", {
  gsd <- gsDesign(
    test.type = 4, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
    delta = delta, alpha = alpha, beta = beta)
  gsdv <- gs_design_npe(
    theta = delta, info = timing, beta = beta,
    binding = FALSE, # Use this for test.type=4 and 6
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_spending_bound,
    lpar = list(sf = sfl, total_spend = beta, param = sflpar)
    )%>% 
    filter(hypothesis == "H1")

  expect_equal(gsd$upper$bound, gsdv$Z[1:K], tolerance = 7e-6)
  expect_equal(gsd$lower$bound, gsdv$Z[(K+1):(2*K)], tolerance = 9e-6)
  expect_equal(gsd$n.I, gsdv$info[1:K], tolerance = .04) # While tolerance should not be problematic, it seems large
  
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(
    theta = 0, info = gsdv$info[1:K],
    upper = gs_spending_bound,
    upar = list(sf = sfu, total_spend = alpha, param = sfupar),
    lower = gs_b,
    lpar = rep(-Inf, K)
    ) %>% 
    filter(Bound == "Upper" & hypothesis == "H1")
  expect_equal(gsdv0$Probability[1:K], sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})


test_that("Two-sided asymmetric design fails to reproduce gsDesign test.type=5 bounds", {
  astar <- 0.2
  gsd <- gsDesign(test.type = 5, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
                  delta = delta, alpha = alpha, beta = beta, astar = astar)
  gsdv <- gs_design_npe(theta = delta, info = timing, beta = beta,
                        theta1 = 0, # Spending for lower bound under H0
                        binding = TRUE, # Use this for test.type=3 and 5
                        upper = gs_spending_bound,
                        upar = list(sf = sfu, total_spend = alpha, param = sfupar),
                        lower = gs_spending_bound,
                        lpar = list(sf = sfl, total_spend = astar, param = sflpar)
  )

  expect_equal(gsd$upper$bound, gsdv$Z[1:K], tolerance = 7e-6)
  expect_equal(gsd$lower$bound, gsdv$Z[(K+1):(2*K)], tolerance = 9e-6)
  expect_equal(gsd$n.I, gsdv$info[1:K], tolerance = .04) # While tolerance should not be problematic, it seems large
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(theta = 0, info = gsdv$info[1:K],
                        upper = gs_spending_bound,
                        upar = list(sf = sfu, total_spend = alpha, param = sfupar),
                        lower = gs_spending_bound,
                        lpar = list(sf = sfu, total_spend = alpha, param = sfupar)
  )
  expect_equal(gsdv0$Probability[1:K], sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})

test_that("Two-sided asymmetric design fails to reproduce gsDesign test.type=6 bounds", {
  astar <- 0.2
  gsd <- gsDesign(test.type = 6, k = K, sfu = sfu, sfupar = sfupar, sfl = sfl, sflpar = sflpar, timing = timing,
                  delta = delta, alpha = alpha, beta = beta, astar = astar)
  gsdv <- gs_design_npe(theta = delta, info = timing, beta = beta,
                        theta1 = 0, # Spending for lower bound under H0
                        binding = FALSE, # Use this for test.type=3 and 5
                        upper = gs_spending_bound,
                        upar = list(sf = sfu, total_spend = alpha, param = sfupar),
                        lower = gs_spending_bound,
                        lpar = list(sf = sfl, total_spend = astar, param = sflpar)
  )

  expect_equal(gsd$upper$bound, gsdv$Z[1:K], tolerance = 7e-6)
  expect_equal(gsd$lower$bound, gsdv$Z[(K+1):(2*K)], tolerance = 9e-6)
  expect_equal(gsd$n.I, gsdv$info[1:K], tolerance = .04) # While tolerance should not be problematic, it seems large
  # get design properties under null hypothesis (theta = 0)
  gsdv0 <- gs_power_npe(theta = 0, info = gsdv$info[1:K],
                        upper = gs_spending_bound,
                        upar = list(sf = sfu, total_spend = alpha, param = sfupar),
                        lower = gs_spending_bound,
                        lpar = list(sf = sfu, total_spend = alpha, param = sfupar)
  )
  expect_equal(gsdv0$Probability[1:K], sfu(alpha = alpha, t = timing, param = sfupar)$spend)
})


