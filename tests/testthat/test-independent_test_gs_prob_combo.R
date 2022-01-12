####### 1 analysis scenario#####
lower <- -0.6
upper <- 0.4
rho <- c(1,1,0,0)
gamma <- c(0,1,0,1)
tau <- c(-1,-1,-1,-1)
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 30),
                              rate = c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate=log(2)/c(9, 18),
                            hr = c(.9, .6),
                            dropoutRate = rep(.001, 2))
arm <- gs_create_arm(enrollRates = enrollRates,
                     failRates = failRates,
                     ratio = 1,
                     total_time = 1e6)
sigma <- gs_sigma2_combo(arm0 = arm$arm0,
                         arm1 = arm$arm1,
                         tmax = 30,
                         rho = rho,
                         gamma=gamma,
                         tau = rep(-1, length(rho)),
                         approx="asymptotic")
corr <- cov2cor(sigma)
n_test <- length(rho)
theta <- rep(0, n_test)
analysis <- 1
prob <- gs_prob_combo(lower_bound = rep(lower, n_test),
                      upper_bound = rep(upper, n_test),
                      analysis = analysis,
                      theta = theta,
                      corr = corr,
                      algorithm = GenzBretz(maxpts = 1e5, abseps = 1e-5))

p_efficacy <- pmvnorm_combo(lower = rep(upper, n_test),
                            upper = rep(Inf, n_test),
                            group = analysis,
                            mean = theta,
                            corr = corr)
p_futility <- pmvnorm_combo(lower = rep(-Inf, n_test),
                            upper = rep(lower, n_test),
                            group = analysis,
                            mean = theta,
                            corr = corr)
test_that ( "p efficacy",{
  expect_equal(prob$Probability[1], p_efficacy[1], tolerance=0.001)
})
test_that ( "p futility",{
  expect_equal(prob$Probability[2], p_futility[1], tolerance=0.001)
})

####### 2 analysis scenario#####
lower <- c(-0.2, -0.3)
upper <- c(0.3, 0.4)
rho <- c(1, 1, 0, 0)
gamma <- c(0, 1, 0, 1)
tau <- c(-1, -1, -1, -1)
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 30),
                              rate = c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate = log(2)/c(9, 18),
                            hr = c(.9, .6),
                            dropoutRate = rep(.001, 2))
arm <- gs_create_arm(enrollRates = enrollRates,
                     failRates = failRates,
                     ratio = 1,
                     total_time = 1e6)
sigma <- gs_sigma2_combo (arm0 = arm$arm0,
                          arm1 = arm$arm1,
                          tmax = 30,
                          rho = rho,
                          gamma = gamma,
                          tau = rep(-1, length(rho)),
                          approx ="asymptotic")
corr <- cov2cor(sigma)
n_test <- length(rho)
theta <- rep(0, n_test)
analysis <- c(1,2)
prob <- gs_prob_combo(lower_bound = rep(lower, n_test),
                      upper_bound = rep(upper, n_test),
                      analysis = analysis,
                      theta = theta,
                      corr = corr,
                      algorithm = GenzBretz(maxpts= 1e5, abseps= 1e-5))
c <- c(1,3)
corr1 = corr[c,c]
p_efficacy_1 <- pmvnorm_combo(lower = rep(upper[1], 2),
                              upper = rep(Inf, 2),
                              group = 1,
                              mean = theta[c],
                              corr = corr1)
p_futility_1 <- pmvnorm_combo(lower = rep(-Inf, 2),
                              upper = rep(lower[1], 2),
                              group = 1,
                              mean = theta[c],
                              corr = corr1)
p_efficacy_2 <- pmvnorm_combo(lower=c(lower[1], upper[2]),
                              upper=c(upper[1], Inf),
                              group = analysis,
                              mean = theta,
                              corr = corr)
p_futility_2 <- pmvnorm_combo(lower=c(lower[1], -Inf),
                              upper=c(upper[1], lower[2]),
                              group = analysis,
                              mean = theta,
                              corr = corr)
test_that("p efficacy1",{
  expect_equal(prob$Probability[1], p_efficacy_1[1], tolerance = 0.001)
})
test_that("p futility1",{
  expect_equal(prob$Probability[3], p_futility_1[1], tolerance = 0.001)
})
test_that("p efficacy2",{
  expect_equal(prob$Probability[2], p_efficacy_1[1]+p_efficacy_2[1], tolerance = 0.001)
})
test_that("p futility2",{
  expect_equal(prob$Probability[4], p_futility_1[1]+p_futility_2[1], tolerance = 0.001)
})
