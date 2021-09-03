#--------------------------------------
# Simulation for gsDesign with WLR
#--------------------------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Set up Simulation Environment

# task_id <- 1
set.seed(task_id)

library(gsDesign2)
library(simtrial)
library(dplyr)
library(survival)
library(mvtnorm)
library(survMisc)
library(purrr)
#--------------------------------------
# tenFHcorr from simtrial MagirrBurman branch
#--------------------------------------
tenFHcorr <- function(x=simPWSurv(n=200) %>% cutDataAtCount(100) %>%
                        tensurv(txval = "Experimental"),
                      rg=tibble(rho=c(0,0,1,1),gamma=c(0,1,0,1)),
                      corr=TRUE)
{
  nr <- nrow(rg)
  rhoave <- (matrix(rg$rho,nrow=nr,ncol=nr)+matrix(rg$rho,nrow=nr,ncol=nr,byrow=TRUE))/2
  gamave <- (matrix(rg$gamma,nrow=nr,ncol=nr)+matrix(rg$gamma,nrow=nr,ncol=nr,byrow=TRUE))/2
  # Convert back to tibble
  rg2 <- tibble(rho=as.numeric(rhoave), gamma=as.numeric(gamave))
  rgu <- rg2 %>% unique()
  rgFH <- rg2 %>% left_join(tenFH(x,rgu,returnVariance=TRUE),by=c("rho"="rho","gamma"="gamma"))
  Z <- rgFH$Z[(0:(nrow(rg)-1))*nrow(rg)+1:nrow(rg)]
  c <- matrix(rgFH$Var,nrow=nrow(rg),byrow=TRUE)
  if (corr) c <- stats::cov2cor(c)
  names(c) <- paste("V",1:ncol(c),sep="")
  cbind(rg,Z,as_tibble(c))
}


sim_gsd_pMaxCombo <- function(N = ceiling(454.60),
                              enrollRates = tibble::tibble(Stratum = "All", duration = 12, rate = N/12),
                              failRates = tibble::tibble(Stratum = "All",
                                                         duration = c(4, 100),
                                                         failRate = log(2) / 15,  # median survival 15 month
                                                         hr = c(1, .6),
                                                         dropoutRate = 0.001),
                              fh_test,
                              # lower = c(-0.8644507,  1.2780679,  2.1082168), #for example 2
                              # upper = c(3.822543, 2.628201, 2.110026)
                              lower, #for example 1
                              upper

){

  # analysis time
  time <- sort(unique(fh_test$analysisTimes))

  # all tests
  u_fh_test <- unique(fh_test[, c("rho", "gamma")])

  failRates0 <- tibble::tibble(Stratum    = failRates$Stratum,
                               period     = 1:nrow(failRates),
                               Treatment  = "Control",
                               duration   = failRates$duration,
                               rate       = failRates$failRate)

  failRates1 <- tibble::tibble(Stratum    = failRates$Stratum,
                               period     = 1:nrow(failRates),
                               Treatment  = "Experimental",
                               duration   = failRates$duration,
                               rate       = failRates$failRate * failRates$hr)

  dropoutRates0 = tibble::tibble(Stratum = failRates$Stratum,
                                 period  = 1:nrow(failRates),
                                 Treatment = "Control",
                                 duration = failRates$duration,
                                 rate = failRates$dropoutRate)

  dropoutRates1 = tibble::tibble(Stratum = failRates$Stratum,
                                 period  = 1:nrow(failRates),
                                 Treatment = "Experimental",
                                 duration = failRates$duration,
                                 rate = failRates$dropoutRate)



  sim <- simtrial::simPWSurv(n = as.numeric(N),
                             enrollRates  = enrollRates,
                             failRates    = bind_rows(failRates0, failRates1),
                             dropoutRates = bind_rows(dropoutRates0, dropoutRates1))

  # Analysis for each interim analysis
  foo <- function(t,sim){
    sim_cut <- sim %>% simtrial::cutData(cutDate = t)

    # Total events
    d <- sum(sim_cut$event)

    # Weighted log rank test
    # z0 <- c()
    # for(i in 1:nrow(u_fh_test)){
    #   z0[[i]] <- sim_cut %>% tensurv(txval = "Experimental") %>% tenFH(rg = u_fh_test[i, ])
    # }
    # z0 <- bind_rows(z0)

    z <- sim_cut %>% tensurv(txval = "Experimental") %>% tenFHcorr(rg = u_fh_test)
    #pMC = pMaxCombo(z)
    bind_cols(n = N,
              analysisTimes = t,
              d = d,
              z = z[, c("rho", "gamma", "Z")], # pMaxCombo = pMC
    )
  }

  res <- bind_rows(lapply(time, foo, sim = sim)) # for example 2

  fh_res <- merge(fh_test, res, all = TRUE)


  # sequential test procedure
  fh_res <- fh_res %>% subset(! is.na(Analysis)) %>%
                        group_by(n, d, analysisTimes, Analysis) %>%
                        summarise(z = max( - Z))
  z <- fh_res$z

  p <- ! (z < upper & z > lower)

  test_lower <- rep(FALSE, length(time))
  test_upper <- rep(FALSE, length(time))

  test_i <- which(p)[1]
  if(z[test_i] > upper[test_i]){
    test_upper[test_i] <- TRUE
  }

  if(z[test_i] < lower[test_i]){
    test_lower[test_i] <- TRUE
  }

  fh_res$lower <- test_lower
  fh_res$upper <- test_upper
  fh_res$lower_bound <- lower
  fh_res$upper_bound <- upper

  list(test = res, res = fh_res)
}


fh_test1 <- rbind( data.frame(rho = 0, gamma = 0, tau = -1,
                             test = 1,
                             Analysis = 1:3,
                             analysisTimes = c(12, 24, 36)),
                  data.frame(rho = c(0, 0.5), gamma = 0.5, tau = -1,
                             test = 2:3,
                             Analysis = 3, analysisTimes = 36)
)

fh_test2 <- data.frame(rho = c(0, 0), gamma = c(0, 0.5), tau = -1,
                       analysisTimes = rep(c(12, 24, 36), each = 2),
                       Analysis = rep(1:3, each = 2),
                       test = rep(1:2, 3))
results <- list(
  s01 = sim_gsd_pMaxCombo(N = ceiling(444.78),
                          upper = c(3.71, 2.51, 1.99),
                          lower = c(-0.24, 1.17, 1.99),
                          fh_test = fh_test1),

  s02 = sim_gsd_pMaxCombo(N = ceiling(348.2),
                          upper = c(3.71, 2.51, 1.99),
                          lower = c(-0.24, 1.17, 1.99),
                          fh_test = fh_test2),

  s03 = sim_gsd_pMaxCombo(N = ceiling(301.15),
                          upper = c(6.18, 2.80, 2.10),
                          lower = c(-2.72, 0.65, 2.10),
                          fh_test = fh_test1)
)

result_res  <- ungroup(bind_rows(map(results, "res"), .id = "scenario"))
result_test <- ungroup(bind_rows(map(results, "test"), .id = "scenario"))

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(result_res, result_test, file = filename)

#----------------------
# HPC Submission code
#----------------------

# cd /SFS/scratch/zhanyilo/simu_gsd_combo
# rm *
# module add R/4.0.2
# qsub -t 1:10000 ~/runr.sh ~/gsdmvn/simulation/simu_gsd_combo.R

#-------------------------------
# Summarize simulation results
#-------------------------------
library(purrr)
library(dplyr)
library(tidyr)

path <- "/SFS/scratch/zhanyilo/simu_gsd_combo/"

res <- list()
test <- list()
for(i in 1:10000){
  load(file.path(path, paste0(i, ".Rdata")))
  try({
    res[[i]]  <- result_res
    test[[i]] <- result_test
  })
}


res <- bind_rows(res)
test <- bind_rows(test)


sim_res_combo <-
  res %>% group_by(scenario, n, Analysis) %>%
        summarise_all(mean) %>%
        group_by(scenario, n) %>%
        mutate(lower = cumsum(lower),
               upper = cumsum(upper)) %>% data.frame()

library(tidyr)
corr1 <- test %>% subset(scenario == "s01") %>%
                  as_tibble() %>%
                  select(-scenario, -d, -n) %>%
                  unite("name", rho, gamma, analysisTimes) %>%
                  group_by(name) %>%
                  mutate(id = 1:n()) %>%
                  ungroup() %>%
                  arrange(name) %>%
                  pivot_wider(names_from = name,
                              values_from = Z) %>%
                  select(- id) %>%
                  cor()

corr2 <- test %>% subset(scenario == "s02") %>%
                  as_tibble() %>%
                  select(-scenario, -d, -n) %>%
                  unite("name", rho, gamma, analysisTimes) %>%
                  group_by(name) %>%
                  mutate(id = 1:n()) %>%
                  ungroup() %>%
                  arrange(name) %>%
                  pivot_wider(names_from = name,
                              values_from = Z) %>%
                  select(- id) %>%
                  cor()

save(sim_res_combo, corr1, corr2, file = "./simulation/simu_gsd_combo.Rdata")
