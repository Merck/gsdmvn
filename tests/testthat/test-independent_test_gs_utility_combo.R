#####log-rank multiple analysis#####
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 30),
                              rate = c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate = log(2)/c(9, 18),
                            hr = c(.9, .6),
                            dropoutRate = rep(.001, 2))
analysisTimes <- c(12, 24, 36)
n_analysis <- length(analysisTimes)
fh_test <- rbind(data.frame(rho = 0,
                            gamma = 0,
                            tau = -1,
                            test = 1,
                            Analysis = 1:3,
                            analysisTimes = analysisTimes)
)
gs_arm <- gsdmvn:::gs_create_arm(enrollRates,
                                 failRates,
                                 ratio = 1,                       # Randomization ratio
                                 total_time = max(analysisTimes)) # Total study duration

utility_combo <- gsdmvn:::gs_utility_combo(enrollRates = enrollRates,
                                           failRates = failRates,
                                           fh_test = fh_test,
                                           ratio = 1,
                                           algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5))

info_combo_test <- gsdmvn:::gs_info_combo(enrollRates = enrollRates,
                                          failRates = failRates,
                                          ratio = 1,
                                          analysisTimes = analysisTimes,
                                          rho = 0,
                                          gamma = 0)
test_that("gs_utility_combo output correct info as gs_info_combo",{
  expect_equal(utility_combo$info[1:11], info_combo_test[1:11])
})

theta_test<- (- info_combo_test$delta) / sqrt(info_combo_test$sigma2)
test_that("gs_utility_combo output correct theta effect as gs_info_combo",{
  expect_equal(utility_combo$theta, theta_test)
})

info <- info_combo_test[[10]]
cov <- matrix(0, n_analysis, n_analysis)
for (i in 1:n_analysis) {
  for (j in 1:n_analysis) {
    k <- min(i,j)
    cov[i,j] <- info[k]/(info[i]*info[j])
  }
}
corr_test <- cov2cor(cov)
test_that ( "gs_utility_combo output correct correlation matrix as gs_info_combo",{
  expect_equal(utility_combo$corr, corr_test)
})


##### multiple test analysis#####
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 30),
                              rate = c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate = log(2)/c(9, 18),
                            hr = c(.9, .6),
                            dropoutRate = rep(.001, 2))
analysisTimes <- 36
n_analysis <- length(analysisTimes)
rho <- c(0, 0.5, 1)
gamma <- c(0.5, 0.5, 0.5)
tau <- c(-1, -1, -1)
fh_test <- rbind(data.frame(rho = rho,
                            gamma = gamma,
                            tau = tau,
                            test = 1:3,
                            Analysis = 1,
                            analysisTimes = analysisTimes)
)
gs_arm <- gsdmvn:::gs_create_arm(enrollRates,
                                 failRates,
                                 ratio = 1,                   # Randomization ratio
                                 total_time = max(analysisTimes)) # Total study duration

utility_combo <- gsdmvn:::gs_utility_combo(enrollRates = enrollRates,
                                           failRates = failRates,
                                           fh_test = fh_test,
                                           ratio = 1,
                                           algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5))

info_combo_test <- gsdmvn:::gs_info_combo(enrollRates = enrollRates,
                                          failRates = failRates,
                                          ratio = 1,
                                          analysisTimes = analysisTimes,
                                          rho = rho,
                                          gamma = gamma)
test_that ( "gs_utility_combo output correct info as gs_info_combo",{
  expect_equal(utility_combo$info[1:11], info_combo_test[1:11])
})

theta_test<- (- info_combo_test$delta) / sqrt(info_combo_test$sigma2)
test_that ( "gs_utility_combo output correct theta effect as gs_info_combo",{
  expect_equal(utility_combo$theta, theta_test)
})

sigma2 <- gsdmvn:::gs_sigma2_combo(arm0 = gs_arm$arm0,
                                   arm1 = gs_arm$arm1,
                                   tmax = analysisTimes,
                                   rho = rho,
                                   gamma = gamma,
                                   tau = tau)
corr_test <-cov2cor(sigma2)
test_that ( "gs_utility_combo output correct correlation matrix as gs_info_combo",{
  expect_equal(utility_combo$corr, corr_test)
})
