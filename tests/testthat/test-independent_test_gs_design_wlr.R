####NOTE: all reference numbers come from simulation results in
###https://keaven.github.io/gsd-deming/wlr.html#wlr
###please look for sections titled "Simulation results based on 10,000 replications."

test_that("Validate the function based on examples with simulation results",{
  x <- gsDesign::gsSurv(
    k = 3, test.type = 4, alpha = 0.025,
    beta = 0.1, astar = 0, timing = c(1),
    sfu = gsDesign::sfLDOF, sfupar = c(0),
    sfl = gsDesign::sfLDOF, sflpar = c(0),
    lambdaC = c(0.1),
    hr = 0.6, hr0 = 1, eta = 0.01,
    gamma = c(10),
    R = c(12), S = NULL,
    T = 36, minfup = 24, ratio = 1
  )
  enrollRates <- tibble::tibble(Stratum = "All",
                                duration = 12,
                                rate = 500 / 12)
  failRates <- tibble::tibble(Stratum = "All",
                              duration = c(4, 100),
                              failRate = log(2) / 15, # Median survival 15 month
                              hr = c(1, 0.6),
                              dropoutRate = 0.001
  )
  ## Randomization Ratio is 1:1
  ratio <- 1
  ## Type I error (one-sided)
  alpha <- 0.025
  ## Power (1 - beta)
  beta <- 0.2
  power <- 1 - beta
  # Interim Analysis Time
  analysisTimes <- c(12, 24, 36)
  #logrank test
  lrk <- gsdmvn::gs_design_wlr(
    enrollRates = enrollRates,
    failRates = failRates,
    weight = function(x, arm0, arm1) {
      gsdmvn::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0)
    },
    ratio = ratio, alpha = alpha, beta = beta,
    upar = x$upper$bound,
    lpar = x$lower$bound,
    analysisTimes = c(12, 24, 36)
  )$bounds %>%
    dplyr::mutate_if(is.numeric, round, digits = 2) %>%
    select(Analysis, Bound, Time, N, Events, AHR, Probability) %>%
    tidyr::pivot_wider(names_from = Bound, values_from = Probability)
  #FH(0,1)
  fh01 <- gsdmvn::gs_design_wlr(
    enrollRates = enrollRates,
    failRates = failRates,
    weight = function(x, arm0, arm1) {
      gsdmvn::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 1)
    },
    ratio = ratio, alpha = alpha, beta = beta,
    upar = x$upper$bound,
    lpar = x$lower$bound,
    analysisTimes = analysisTimes
  )$bounds %>%
    dplyr::mutate_if(is.numeric, round, digits = 2) %>%
    select(Analysis, Bound, Time, N, Events, AHR, Probability) %>%
    tidyr::pivot_wider(names_from = Bound, values_from = Probability)
  #FH(0,0.5)
  fh0d5 <- gsdmvn::gs_design_wlr(
    enrollRates = enrollRates,
    failRates = failRates,
    weight = function(x, arm0, arm1) {
      gsdmvn::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)
    },
    ratio = ratio, alpha = alpha, beta = beta,
    upar = x$upper$bound,
    lpar = x$lower$bound,
    analysisTimes = analysisTimes
  )$bounds %>%
    dplyr::mutate_if(is.numeric, round, digits = 2) %>%
    select(Analysis, Bound, Time, N, Events, AHR, Probability) %>%
    tidyr::pivot_wider(names_from = Bound, values_from = Probability)
  #FH(0.5,0.5)
  fh5d5 <- gsdmvn::gs_design_wlr(
    enrollRates = enrollRates,
    failRates = failRates,
    weight = function(x, arm0, arm1) {
      gsdmvn::wlr_weight_fh(x, arm0, arm1, rho = 0.5, gamma = 0.5)
    },
    ratio = ratio, alpha = alpha, beta = beta,
    upar = x$upper$bound,
    lpar = x$lower$bound,
    analysisTimes = analysisTimes
  )$bounds %>%
    dplyr::mutate_if(is.numeric, round, digits = 2) %>%
    select(Analysis, Bound, Time, N, Events, AHR, Probability) %>%
    tidyr::pivot_wider(names_from = Bound, values_from = Probability)

  #logrank part
  expect_equal(object = as.numeric(lrk$N), expected = rep(386, 3), tolerance = 3)
  expect_equal(object = as.numeric(lrk$Events), expected = c(82.77, 190.05, 255.61), tolerance = 3)
  expect_equal(object = as.numeric(lrk$AHR), expected = c(0.87, 0.72, 0.69), tolerance = 0.3)
  expect_equal(object = as.numeric(lrk$Upper), expected = c(0.00, 0.41, 0.80), tolerance = 0.3)
  expect_equal(object = as.numeric(lrk$Lower), expected = c(0.07, 0.14, 0.20), tolerance = 0.3)
  #fh(0,1)
  expect_equal(object = as.numeric(fh01$N), expected = rep(317, 3), tolerance=3)
  expect_equal(object = as.numeric(fh01$Events), expected = c(68.01, 156.13, 210.06), tolerance = 3)
  expect_equal(object = as.numeric(fh01$AHR), expected = c(0.76, 0.65, 0.63), tolerance = 0.3)
  expect_equal(object = as.numeric(fh01$Upper), expected = c(0.00, 0.45, 0.79), tolerance = 0.3)
  expect_equal(object = as.numeric(fh01$Lower), expected = c(0.04, 0.12, 0.21), tolerance = 0.3)
  #fh(0,0.5)
  expect_equal(object = as.numeric(fh0d5$N), expected = rep(314, 3), tolerance = 3)
  expect_equal(object = as.numeric(fh0d5$Events), expected = c(67.21, 154.43, 207.92), tolerance = 3)
  expect_equal(object = as.numeric(fh0d5$AHR), expected = c(0.81, 0.67, 0.65), tolerance = 0.3)
  expect_equal(object = as.numeric(fh0d5$Upper), expected = c(0.00, 0.44, 0.79), tolerance = 0.3)
  expect_equal(object = as.numeric(fh0d5$Lower), expected = c(0.05, 0.12, 0.21), tolerance = 0.3)
  #fh(0.5,0.5)
  expect_equal(object = as.numeric(fh5d5$N), expected = rep(317, 3), tolerance = 3)
  expect_equal(object = as.numeric(fh5d5$Events), expected = c(67.87, 155.86, 209.82), tolerance = 3)
  expect_equal(object = as.numeric(fh5d5$AHR), expected = c(0.81, 0.68, 0.66), tolerance = 0.3)
  expect_equal(object = as.numeric(fh5d5$Upper), expected = c(0.00, 0.43, 0.80), tolerance = 0.3)
  expect_equal(object = as.numeric(fh5d5$Lower), expected = c(0.06, 0.12, 0.20), tolerance = 0.3)
})
