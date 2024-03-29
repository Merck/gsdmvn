# Test 1: compare results with gsDesign2::AHR ####

testthat::test_that("compare results with gsDesign2::AHR in the situation of single analysis",{
  enrollRates <- tibble::tibble(Stratum = "All",
                                duration = c(2, 2, 10),
                                rate = c(3, 6, 9))
  failRates <- tibble::tibble(Stratum = "All",
                              duration = c(3, 100),
                              failRate = log(2)/c(9, 18),
                              hr = c(0.9, 0.6),
                              dropoutRate = rep(0.001, 2))
  totalDuration <- 30
  analysisTimes <- totalDuration

  out <- gs_design_ahr(enrollRates = enrollRates,
                       failRates = failRates,
                       analysisTimes = analysisTimes)

  testthat::expect_equal(out$bounds %>% select(Time, AHR),
                         gsDesign2::AHR(enrollRates = enrollRates,
                                        failRates = failRates,
                                        totalDuration = totalDuration) %>% select(Time, AHR))

  #update enrollRates for gsDesign2::AHR to make Events/info/info0 also match in outputs
  enrollRates1 <- enrollRates %>% mutate(rate = rate * c(out$bounds$N / (duration %*% rate)))

  testthat::expect_equal(out$bounds %>% select(Time, AHR, Events, info, info0),
                         gsDesign2::AHR(enrollRates = enrollRates1,
                                        failRates = failRates,
                                        totalDuration = totalDuration) %>% select(Time, AHR, Events, info, info0))
  })

testthat::test_that("compare results with gsDesign2::AHR in the situation with IF and multiple analysis times specified",{
  enrollRates <- tibble::tibble(Stratum = "All",
                                duration = c(2, 2, 10),
                                rate = c(3, 6, 9))
  failRates <- tibble::tibble(Stratum = "All",
                              duration = c(3, 100),
                              failRate = log(2)/c(9, 18),
                              hr = c(0.9, 0.6),
                              dropoutRate = rep(0.001, 2))
  totalDuration <- c(12, 25, 36)
  analysisTimes <- totalDuration

  out <- gs_design_ahr(enrollRates = enrollRates,
                       failRates = failRates,
                       analysisTimes = analysisTimes)

  testthat::expect_equal(out$bounds %>% select(Time, AHR) %>% dplyr::distinct (.keep_all = TRUE),
                         gsDesign2::AHR(enrollRates = enrollRates,
                                        failRates = failRates,
                                        totalDuration = totalDuration) %>% select(Time, AHR))

  #update enrollRates for gsDesign2::AHR to make Events/info/info0 also match in outputs
  enrollRates1 <- enrollRates %>% mutate(rate = rate * c(max(out$bounds$N) / (duration %*% rate)))

  testthat::expect_equal(out$bounds %>%
                           select(Time, AHR, Events, info, info0) %>%
                           dplyr::distinct (.keep_all = TRUE),
                         gsDesign2::AHR(enrollRates = enrollRates1,
                                        failRates = failRates,
                                        totalDuration = totalDuration) %>%
                           select(Time, AHR, Events, info, info0))
})
