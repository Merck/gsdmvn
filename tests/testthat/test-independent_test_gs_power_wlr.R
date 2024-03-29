test_that("Check using gs_info_wlr and gs_power_npe", {
  enrollRates <- tibble::tibble(Stratum = "All",
                                duration = 12,
                                rate = 500 / 12)
  failRates <- tibble::tibble(
    Stratum = "All",
    duration = c(4, 100),
    failRate = log(2) / 15, # Median survival 15 months
    hr = c(1, .6), # Delay effect after 4 months
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

  #create arms
  # Define study design object in each arm
  gs_arm <- gsdmvn:::gs_create_arm(
    enrollRates,
    failRates,
    ratio = 2, # Randomization ratio
    total_time = 36 # Total study duration
  )
  arm0 <- gs_arm[["arm0"]]
  arm1 <- gs_arm[["arm1"]]
  #calculate all pieces of information
  weight <- function(x, arm0, arm1) {
    gsdmvn::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 1)
  }
  gs_info <- gsdmvn::gs_info_wlr(
    enrollRates,
    failRates,
    ratio,
    analysisTimes = analysisTimes,
    weight = weight
  )
  fh01 <-gs_info %>% dplyr::mutate_if(is.numeric, round, digits = 5)

  up <-gsDesign(k = length(fh01$Events),
                test.type = 1,
                n.I = fh01$Events,
                maxn.IPlan = max(fh01$Events),
                sfu = sfLDOF,
                sfupar = NULL)$upper$bound

  npe <- gs_power_npe(theta = fh01$theta,
                      info = fh01$info,
                      info0 = fh01$info0,
                      binding = F,
                      upper = gs_b,
                      lower = gs_b,
                      upar = up,
                      lpar = c(qnorm(.1), rep(-Inf, length(fh01$Events) - 1)),
                      test_upper = T,
                      test_lower = T,
                      r = 18,
                      tol = 1e-6)

 #output
  gspow <- gs_power_wlr(enrollRates = enrollRates,
                        failRates = failRates,
                        ratio = ratio,               # Experimental:Control randomization ratio
                        weight = weight,
                        approx = "asymptotic",
                        events = fh01$Events, # Targeted events of analysis
                        analysisTimes = NULL,   # Targeted times of analysis
                        binding = FALSE,
                        upper = gs_b, # Default is Lan-DeMets approximation of
                        upar = up,
                        lower = gs_b,
                        lpar = c(qnorm(.1), rep(-Inf, length(fh01$Events) - 1)), # Futility only at IA1
                        test_upper = TRUE,
                        test_lower = TRUE,
                        r = 18,
                        tol = 1e-6)

  #tests
  expect_equal(object = as.numeric(gspow$Time), expected = rep(fh01$Time, 2), tolerance = 0.0001)
  expect_equal(object = as.numeric(gspow$Events), expected = rep(fh01$Events, 2), tolerance = 1)
  expect_equal(object = as.numeric(gspow$Z), expected = npe$Z, tolerance = 0.1)
  expect_equal(object = as.numeric(gspow$Probability), expected = npe$Probability, tolerance = 0.001)
  expect_equal(object = as.numeric(gspow$AHR), expected = rep(fh01$AHR, 2), tolerance = 0.01)
  expect_equal(object = as.numeric(gspow$theta), expected = rep(fh01$theta, 2), tolerance = 0.001)
  expect_equal(object = as.numeric(gspow$info), expected = rep(fh01$info, 2), tolerance = 0.001)
  expect_equal(object = as.numeric(gspow$info0), expected = rep(fh01$info0, 2), tolerance = 0.001)
})
