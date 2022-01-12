
# \\\\test get_combo_weight\\\\\
rho <- c(1, 1, 0, 0)
gamma <- c(0, 1, 0, 1)
tau <- c(-1, -1, -1, -1)

weight <- get_combo_weight(rho, gamma, tau)
weight1_rho <- substring(weight[[1]], 125, 130)
weight2_rho <- substring(weight[[2]], 125, 130)
weight1_gamma <- substring(weight[[1]], 133, 140)
weight2_gamma <- substring(weight[[2]], 133, 140)
weight1_tau <- substring(weight[[1]], 143, 148)
weight2_tau <- substring(weight[[2]], 143, 148)

test_that ( "get_combo_weight output correct rho1",{
  expect_equal(weight1_rho, "rho =1")
})
test_that ( "get_combo_weight output correct rho2",{
  expect_equal(weight2_rho, "rho =1")
})

test_that ( "get_combo_weight output correct gamma1",{
  expect_equal(weight1_gamma, "gamma =0")
})
test_that ( "get_combo_weight output correct gamma2",{
  expect_equal(weight2_gamma, "gamma =1")
})

# \\\\test get_combo_weight tau not equal to -1\\\\\
rho <- c(1, 1, 0, 0)
gamma <- c(0, 1, 0, 1)
tau <- c(1, 1, 0, 0)
weight <- get_combo_weight(rho, gamma, tau)
weight1_tau <- substring(weight[[1]], 143, 148)
weight3_tau <- substring(weight[[3]], 143, 148)
test_that ( "get_combo_weight output correct tau1",{
  expect_equal(weight1_tau, "tau =1")
})
test_that ( "get_combo_weight output correct tau3",{
  expect_equal(weight3_tau, "tau =0")
})

# \\\\test gs_delta_combo\\\\\
rho <- c(1,1,0,0)
gamma <- c(0,1,0,1)
tau <- c(-1,-1,-1,-1)
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2, 2, 30),
                              rate = c(3, 6, 9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3, 100),
                            failRate = log(2)/c(9, 18),
                            hr = c(.9,.6),
                            dropoutRate = rep(.001, 2))
arm <- gs_create_arm(enrollRates,
                     failRates,
                     ratio = 1,
                     total_time = 1e6)
delta <- gs_delta_combo(arm0 = arm$arm0,
                        arm1 = arm$arm1,
                        tmax = 30,
                        rho = rho,
                        gamma = gamma,
                        tau = rep(-1, length(rho)),
                        approx = "asymptotic",
                        normalization = FALSE)
for (i in 1:4) {
  weight_test1 <- get_combo_weight(rho[i], gamma[i], tau[i])
  delta_test1 <- gs_delta_wlr(arm0 = arm$arm0,
                              arm1 = arm$arm1,
                              tmax = 30,
                              weight = eval(parse(text = weight_test1)),
                              approx = "asymptotic", normalization = FALSE)
  test_that ( "gs_delta_combo correctly use gs_delta_wlr 1",{
    expect_identical(delta[i], delta_test1)
  })
}

# \\\\test gs_sigma2_combo\\\\\
sigma2 <- gs_sigma2_combo(arm0 = arm$arm0,
                          arm1 = arm$arm1,
                          tmax = 30,
                          rho = rho,
                          gamma = gamma,
                          tau = rep(-1, length(rho)),
                          approx = "asymptotic")
rho1 <- outer(rho, rho, function(x,y) (x+y)/2 )
gamma1 <- outer(gamma, gamma, function(x,y) (x+y)/2 )
for (i in 1:4) {
  for (j in 1:4) {
    weight_test_ij <- get_combo_weight(rho1[i,j], gamma1[i,j],tau[i])
    sigma_ij=gs_sigma2_wlr(arm0 = arm$arm0,
                           arm1 = arm$arm1,
                           tmax = 30,
                           weight = eval(parse(text = weight_test_ij)) ,
                           approx = "asymptotic")
    test_that ( "gs_sigma2_combo correctly use gs_sigma2_wlr 1",{
      expect_equal(sigma2[i,j], sigma_ij)
    })
  }
}


# \\\\test gs_info_combo\\\\\
rho <- c(1,1,0,0)
gamma <- c(0,1,0,1)
tau <- c(-1,-1,-1,-1)
enrollRates <- tibble::tibble(Stratum = "All",
                              duration = c(2,2,30),
                              rate = c(3,6,9))
failRates <- tibble::tibble(Stratum = "All",
                            duration = c(3,100),
                            failRate = log(2)/c(9,18),
                            hr = c(.9,.6),
                            dropoutRate = rep(.001,2))
info_combo <- gs_info_combo(enrollRates = enrollRates,
                            failRates = failRates,
                            ratio = 1,                # Experimental:Control randomization ratio
                            events = NULL,            # Events at analyses
                            analysisTimes = 30,       # Times of analyses
                            rho = rho,
                            gamma = gamma,
                            tau = rep(-1, length(rho)),
                            approx = "asymptotic"
)
for (i in 1:4) {
  weight_test_i <- get_combo_weight(rho[i], gamma[i], tau[i])
  info_wlr <- gs_info_wlr(enrollRates = enrollRates,
                          failRates = failRates,
                          ratio = 1,                  # Experimental:Control randomization ratio
                          events = NULL,              # Events at analyses
                          analysisTimes = 30,         # Times of analyses
                          weight = eval(parse(text = weight_test_i)),
                          approx = "asymptotic")
  test_that ( "gs_info_combo correctly use gs_info_wlr 1",{
    expect_equal(info_combo$info[i], info_wlr$info[1])
  })
}
