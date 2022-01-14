lower <- -Inf
upper <- 0
mean <- 0.3
n_test <- 4
rho <- c(1, 1, 0, 0)
gamma <- c(0, 1, 0, 1)
tau <- c(-1, -1, -1, -1)

enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 10),
                              rate=c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate = log(2)/c(9, 18),
                            hr = c(.9, .6),
                            dropoutRate = rep(.001, 2))
arm <- gsdmvn:::gs_create_arm(enrollRates = enrollRates,
                              failRates = failRates,
                              ratio = 1,
                              total_time = 1e6)
sigma <- gsdmvn:::gs_sigma2_combo(arm0 = arm$arm0,
                                  arm1 = arm$arm1,
                                  tmax = 30,
                                  rho = rho,
                                  gamma = gamma,
                                  tau = rep(-1, length(rho)),
                                  approx = "asymptotic")
corr <- cov2cor(sigma)

p <- pmvnorm_combo(lower = rep(lower, n_test),
                   upper = rep(upper, n_test),
                   group =2,
                   mean = rep(mean, n_test),
                   corr = corr,
                   algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5))

p_test <- mvtnorm::pmvnorm(lower = rep(lower, n_test),
                           upper = rep(upper, n_test),
                           mean = rep(mean, n_test),
                           corr = corr,
                           sigma = NULL,
                           algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5))

test_that ( "pmvnorm_comb calculate p for One test for all group or lower bound is -Inf.",{
  expect_equal(p[1], p_test[1], tolerance = 0.001)
})

