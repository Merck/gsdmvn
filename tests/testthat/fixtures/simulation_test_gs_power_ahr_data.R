library(gsDesign)
library(simtrial)

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

sampleSize <- ceiling(max(x$eNC) + max(x$eNE))

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

set.seed(1234)
xx <- simfix(nsim=100000,
             sampleSize= sampleSize,
             strata = tibble::tibble(Stratum = "All", p = 1),
             targetEvents = x$n.I[1], # 1st IA only
             block = c(rep("Control", 2), rep("Experimental", 2)),
             enrollRates = tibble::tibble(Stratum = "All",
                                          duration = c(2, 2, 2, 6), rate = c(6, 12, 18, 24)),
             failRates = tibble::tibble(Stratum = "All", duration = 1, failRate =log(2)/9, hr = 0.65,
                                        dropoutRate = 0.001),
             timingType = 2
)

sim.PowerIA <- xx %>% summarize(Power = mean(Z <= -x$upper$bound[1]))

save(sim.PowerIA, file = "./tests/testthat/fixtures/simulation_test_gs_power_ahr_data.Rdata")
