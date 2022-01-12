load("./tests/testthat/fixtures/simulation_test_gs_power_ahr_data.Rdata")

x <- gsSurv(
  k = 2,
  test.type = 1,
  alpha = 0.025,
  beta = 0.2,
  astar = 0,
  timing = 0.7,
  sfu = sfLDOF,
  sfupar = c(0),
  sfl = sfLDOF,
  sflpar = c(0),
  lambdaC = log(2)/9,
  hr = 0.65,
  hr0 = 1,
  eta = 0.001,
  gamma = c(6, 12, 18, 24),
  R = c(2, 2, 2, 6),
  S = NULL,
  T = NULL,
  minfup = NULL,
  ratio = 1
)

#update x with gsDesign() to get integer event counts
x <- gsDesign(
  k = x$k,
  test.type = 1,
  alpha = x$alpha,
  beta = x$beta,
  sfu = x$upper$sf,
  sfupar = x$upper$param,
  n.I = ceiling(x$n.I),
  maxn.IPlan = ceiling(x$n.I[x$k]),
  delta = x$delta,
  delta1 = x$delta1,
  delta0 = x$delta0
)

testthat::test_that("compare power at IA from gs_power_ahr with the simulated result",{
  out <- gs_power_ahr(enrollRates = tibble::tibble(Stratum = "All", duration = c(2, 2, 2, 6), rate = c(6, 12, 18, 24)),
                      failRates = tibble::tibble(Stratum = "All", duration = 1, failRate = log(2)/9, hr = 0.65, dropoutRate = 0.001),
                      ratio = 1,
                      events = x$n.I, #set number of events the same as the design x above from gsDesign()
                      analysisTimes = NULL,
                      binding = FALSE,
                      upper = gs_spending_bound,
                      upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL, theta=0),
                      lower = gs_spending_bound,
                      lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.2, param = NULL, timing = NULL, theta=0),
                      test_upper = TRUE,
                      test_lower = FALSE)

  testthat::expect_equal(out$Probability[1], sim.PowerIA$Power, tolerance = 0.02) #in case test fails, check whether caused by small tolerance
})
