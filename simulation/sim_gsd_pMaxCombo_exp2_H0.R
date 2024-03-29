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


sim_gsd_pMaxCombo_exp2_H0 <- function(N = ceiling(454.60),
                        enrollRates = tibble::tibble(Stratum = "All", duration = 12, rate = N/12),
                        failRates = tibble::tibble(Stratum = "All",
                                                   duration = c(4, 100),
                                                   failRate = log(2) / 15, # median survival 15 month
                                                   hr = c(1, .6),
                                                   dropoutRate = 0.001),
                        rg = tibble(rho = c(0, 0.5), gamma = c(0, 0.5)), # for example 2
                        # rg = list(tibble(rho = c(0), gamma = c(0)), # for example 1
                        #           tibble(rho = c(0), gamma = c(0)),
                        #           tibble(rho = c(0, 0,0.5), gamma = c(0, 0.5, 0.5))),
                        time = c(12, 24, 36),
                        lower = c(-0.8644507,  1.2780679,  2.1082168), #for example 2
                        upper = c(3.822543, 2.628201, 2.110026)
                        #lower = c(-0.2361879,  1.1615317,  2.1203867), #for example 1
                        #upper = c(3.710295, 2.512029, 2.119565)

){


  failRates$failRate <- failRates$failRate * (failRates$hr + 1) / 2
  failRates$hr <- 1
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
    z <- sim_cut %>% tensurv(txval = "Experimental") %>% tenFHcorr(rg = rg)
    #pMC = pMaxCombo(z)
    bind_cols(n = N,
              t = t,
              d = d,
              z = min(z$Z),# pMaxCombo = pMC
              )
  }

  res <- bind_rows(lapply(time, foo, sim = sim)) # for example 2
  # res <- NULL
  # for(t in 1:length(time)){
  #   res = bind_rows(res, foo(time[t], sim, rg = rg[[t]]))
  # }

  names(res) <- tolower(names(res))


  # sequential test procedure
  z <- - res$z
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

  res$lower <- test_lower
  res$upper <- test_upper

  res
}

results <- list(
  s01  = sim_gsd_pMaxCombo_exp2_H0(N = ceiling( 454.60))
)

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(results, file = filename)

#----------------------
# HPC Submission code
#----------------------

# cd /SFS/scratch/shirazam/sim_gsd_pMaxCombo_exp2_H0
# rm *
# cp ~/pkg/gsdmvn/simulation/*.R .
# module add R/4.0.2
# qsub -t 1:10000 ~/runr.sh ~/pkg/gsdmvn/simulation/sim_gsd_pMaxCombo_exp2_H0.R

#-------------------------------
# Summarize simulation results
#-------------------------------

# path <- '/SFS/scratch/shirazam/sim_gsd_pMaxCombo_exp2_H0'
#
#
# res <- list()
# for(i in 1:10000){
#   load(file.path(path, paste0(i, ".Rdata")))
#   try(
#     res[[i]] <- bind_rows(results, .id = "scenario")
#   )
# }
#
# res <- bind_rows(res)
# sim_gsd_pMaxCombo_exp2_H0<-
#   res %>% group_by(scenario, n, t) %>%
#   summarise(events   = mean(d),
#             lower    = mean(lower),
#             upper    = mean(upper)) %>%
#   group_by(scenario, n) %>%
#   mutate(lower = cumsum(lower),
#          upper = cumsum(upper)) %>% data.frame() %>%
#   mutate_if(is.numeric, round, digits = 3)
# save(sim_gsd_pMaxCombo_exp2_H0, file = "./simulation/sim_gsd_pMaxCombo_exp2_H0.Rdata")
