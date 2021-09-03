library(gsdmvn)
library(gsDesign)
library(dplyr)
context("Testing gs_power_npe()")

test_that("Basic call to gs_power_npe fails", {
  # fixed bounds
  x <- try(
  gs_power_npe(theta = .1, info = 1:3, info0 = NULL, binding = FALSE,
       upper=gs_b, lower=gs_b, upar = 4:2, lpar= -1:1,
       test_upper = TRUE, test_lower = TRUE)
  )
  expect_equal(max(class(x)=="data.frame"),1)}
)
test_that("Efficacy bounds for one-sided group sequential spending incorrect", {
# spending bounds
x1s <-
gs_power_npe(theta = .1, info = (1:3)*400, info0 = NULL, binding = FALSE,
             upper= gs_spending_bound,
             upar = list(sf = gsDesign::sfLDOF, param = NULL, total_spend = 0.025),
             lower = gs_b,
             lpar = rep(-Inf, 3)
) %>% filter(Bound == "Upper")
# x1sa <-
#   gs_power_npe(theta = .1, info = (1:3)*400, info0 = NULL, binding = FALSE,
#                upper= gs_spending_bound,
#                upar = list(sf = gsDesign::sfLDOF, param = NULL, total_spend = 0.025),
#                lower = gs_b,
#                lpar = rep(-Inf, 3),
#                tol = 1e-6,
#                r = 31
#   ) %>% filter(Bound == "Upper")
y1s <- gsProbability(k=3, theta= .1, n.I = (1:3) * 400, a = rep(-20, 3),
                     b=gsDesign(k=3, test.type=1, sfu = sfLDOF)$upper$bound)
expect_equal(x1s$Z, y1s$upper$bound, tol = 1e-4)
expect_equal(x1s$Probability, cumsum(y1s$upper$prob), tol = 1e-5)

# gs_power_npe(theta = .1, info = (1:3)*400, info0 = NULL, binding = FALSE,
  #              upper= gs_spending_bound,
  #              upar = list(sf = gsDesign::sfLDOF, param = NULL, total_spend = 0.025),
  #              lower=  gs_spending_bound,
  #              lpar = list(sf = gsDesign::sfLDPocock, param = NULL, total_spend = 0.1),
  #              test_upper = TRUE, test_lower = TRUE)
  # # crossing spending bounds before final analysis (should produce error)
  # gs_power_npe(theta = .1, info = (1:3)*1000, info0 = NULL, binding = FALSE,
  #              upper= gs_spending_bound,
  #              upar = list(sf = gsDesign::sfLDOF, param = NULL, total_spend = 0.025),
  #              lower=  gs_spending_bound,
  #              lpar = list(sf = gsDesign::sfLDPocock, param = NULL, total_spend = 0.1),
  #              test_upper = TRUE, test_lower = TRUE)
  # check vs gsDesign
})
test_that("Comparison of gs_power_npe to gsDesign fails", {
  x <- gsDesign(k = 3, delta = 0.1, test.type = 4, alpha = 0.025, beta = 0.1, timing = c(.5, .75),
                sfu = sfLDOF, sfupar = NULL, sfl = sfHSD, sflpar = -2)
  b1 <- tibble(Analysis = rep(1:3,2), Bound = c(rep("Upper", 3), rep("Lower", 3)),
               Z = c(x$upper$bound,x$lower$bound),
               Probability = c(cumsum(x$upper$prob[,2]), cumsum(x$lower$prob[,2])),
               theta = rep(x$delta, 6),
               info = rep(x$n.I, 2)
  )
  b1a <- gs_power_npe(theta = .1,
                      info = x$n.I,
                      upper= gs_spending_bound,
                      upar = list(sf = gsDesign::sfLDOF, param = NULL, total_spend = 0.025),
                      lower=  gs_spending_bound,
                      lpar = list(sf = gsDesign::sfHSD, param = -2, total_spend = 0.1))
  expect_lt(max(abs(b1$Probability-b1a$Probability)), 2e-6)
})
