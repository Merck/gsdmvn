library(gsdmvn)
library(simtrial)

y <- gs_design_ahr(enrollRates = tibble(Stratum = "All", rate = c( 2.5,5,7.5,10 ), duration = c( 0.125,0.125,0.125,0.125 )),
                  failRates = tibble(Stratum="All", failRate = .2, hr = .5, dropoutRate = .1, duration=1),
                  analysisTimes = 2,
                  upper = gs_b,
                  upar = qnorm(.025),
                  lower = gs_b,
                  lpar = -Inf
)
N <- round(y$bounds$N)

sum.event <- 0
set.seed(1234)
for (i in 1:1000){
  sim <- simPWSurv(n=N,
                   strata=tibble::tibble(Stratum = "All", p = 1),
                   block=c(rep("Control", 2), rep("Experimental", 2)),
                   enrollRates=tibble(Stratum = "All",
                                      rate = c( 2.5,5,7.5,10 )*y$bounds$N/3.125,
                                      duration = c( 0.125,0.125,0.125,0.125 )),
                   failRates=bind_rows(tibble(Stratum="All",period=1,Treatment="Control"     ,duration=1,rate=.2),
                                       tibble(Stratum="All",period=1,Treatment="Experimental",duration=1,rate=.2*.5)),
                   dropoutRates=bind_rows(tibble(Stratum="All" ,period=1,Treatment="Control"     ,duration=100,rate=.1),
                                          tibble(Stratum="All" ,period=1,Treatment="Experimental",duration=100,rate=.1)))

  sum.event <- sum.event + (sim%>%cutData(2)%>%filter(event!=0)%>%group_by(Stratum)%>%count(event))$n
}

sim.Events <- sum.event/1000
save(sim.Events, file = "./tests/testthat/fixtures/simulation_test_gs_design_ahr_data.Rdata")
