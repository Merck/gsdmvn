
testthat::test_that("compare the output of gsdmvn::gs_create_arm with npsurvSS::create_arm()", {

  enrollRates <- tibble::tibble(Stratum = "All",
                                duration = c(2, 2, 10),
                                rate = c(3, 6, 9))

  failRates <- tibble::tibble(Stratum = "All",
                              duration = c(3, 100),
                              failRate = log(2)/c(9, 18),
                              hr = c(.9, .6),
                              dropoutRate = rep(.001, 2))
  ratio <- 1
  totaltime <- 50

  gsarm <- gsdmvn:::gs_create_arm(enrollRates, failRates,
                                 ratio = ratio,
                                 total_time = totaltime)

  expected0 <- npsurvSS::create_arm(size = ratio,
                                    accr_time = sum(enrollRates$duration),
                                    accr_interval = cumsum(enrollRates$duration),
                                    accr_param = enrollRates$rate*enrollRates$duration/sum(enrollRates$rate*enrollRates$duration),
                                    surv_interval = c(0, c(utils::head(failRates$duration, -1), Inf)),
                                    surv_scale = failRates$failRate,
                                    loss_scale =failRates$dropoutRate[1],
                                    total_time = totaltime)



  expected1 <- npsurvSS::create_arm(size = ratio,
                                    accr_time = sum(enrollRates$duration),
                                    accr_interval = cumsum(enrollRates$duration),
                                    accr_param = enrollRates$rate*enrollRates$duration/sum(enrollRates$rate*enrollRates$duration),
                                    surv_interval = c(0, c(utils::head(failRates$duration, -1), Inf)),
                                    surv_scale = failRates$hr * failRates$failRate,
                                    loss_scale = failRates$dropoutRate[1],
                                    total_time = totaltime)

  expect_identical(object = gsarm$arm0, expected = expected0)
  expect_identical(object = gsarm$arm1, expected = expected1)

})
