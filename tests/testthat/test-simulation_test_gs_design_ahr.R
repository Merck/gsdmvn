load("./fixtures/simulation_test_gs_design_ahr_data.Rdata")

testthat::test_that("compare number of events from gs_design_ahr with the simulated result",{
  out <- gs_design_ahr(enrollRates = tibble(Stratum = "All",
                                            rate = c( 2.5, 5, 7.5, 10),
                                            duration = c( 0.125, 0.125, 0.125, 0.125)),
                       failRates = tibble(Stratum = "All",
                                          failRate = .2,
                                          hr = .5,
                                          dropoutRate = .1,
                                          duration=1),
                       analysisTimes = 2,
                       upper = gs_b,
                       upar = qnorm(.025),
                       lower = gs_b,
                       lpar = -Inf
  )

  testthat::expect_equal(out$bounds$Events, sim.Events, tolerance = 0.02) #in case test fails, check whether caused by small tolerance
})
