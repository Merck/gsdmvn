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
library(Matrix)
library(mvtnorm)
library(survMisc)
library(gsdmvn)

#----------------------------------
# Calculate Bound and sample size
#----------------------------------
# enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 500/12)
#
# failRates <- tibble::tibble(Stratum = "All",
#                             duration = c(4, 100),
#                             failRate = log(2) / 15,  # median survival 15 month
#                             hr = c(1, .6),
#                             dropoutRate = 0.001)
#
# ## Randomization Ratio is 1:1
# ratio = 1
#
# ## Type 1 error (one-sided)
# alpha = 0.025
#
# ## Power (1 - beta)
# beta = 0.2
# power = 1 - beta
#
# # Interim Analysis Time
# analysisTimes <- c(12, 24, 36)
#
# rho <- c(0, 0.5, 0, 0)
# gamma <- c(0, 0.5, 0.5, 1)
#
# weight_fun <- list(
#   function(x, arm0, arm1){ gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0)},
#   function(x, arm0, arm1){ gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0.5, gamma = 0.5)},
#   function(x, arm0, arm1){ gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)},
#   function(x, arm0, arm1){ gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 1)}
# )
#
# x_weight <- lapply(weight_fun, function(weight){
#   gsdmvn:::gs_design_wlr(enrollRates = enrollRates, failRates = failRates,
#                          ratio = ratio, alpha = alpha, beta = beta,
#                          weight = weight,
#                          upper = gsdmvn::gs_spending_bound,
#                          lower = gsdmvn::gs_spending_bound,
#                          upar = list(sf = gsDesign::sfLDOF, total_spend = alpha),
#                          lpar = list(sf = gsDesign::sfLDOF, total_spend = beta),
#                          analysisTimes = analysisTimes)$bounds
# })
#
# save(rho, gamma, x_weight, file = "simulation/simu_gsd_wlr_boundary.Rdata")

#--------------------------------------
# maxcombo function for wlr
#--------------------------------------

rm.combo.WLRmax<- function(time        = NULL,
                           status      = NULL,
                           arm         = NULL,
                           wt          = NULL,
                           adjust.methods = c("holm", "hochberg", "hommel", "bonferroni","asymp")[5],
                           ties.method = c("exact", "breslow", "efron")[2],
                           one.sided   = FALSE,
                           HT.est      = FALSE,
                           max         = TRUE,
                           alpha       = 0.025
)
{


  data.anal <- data.frame(time,status,arm)
  fit<- ten(survfit(Surv(time, status) ~ arm, data = data.anal))


  #Testing

  # sink("out1.txt")
  comp(fit, p= sapply(wt, function(x){x[1]}), q= sapply(wt, function(x){x[2]}))
  # sink()

  tst.rslt <- attr(fit, 'lrt')

  #Combination test (exact form)

  if(max & adjust.methods != "asymp"){

    tst.rslt1 <- rbind(tst.rslt[1,],subset(tst.rslt, grepl("FH", tst.rslt$W)))
    Z.tst.rslt1 <- tst.rslt1$Z
    if(one.sided){p.unadjusted <- stats::pnorm(q=tst.rslt1$Z)}
    if(!one.sided){p.unadjusted <- 1- stats::pnorm(q=abs(tst.rslt1$Z)) + stats::pnorm(q=-abs(tst.rslt1$Z))}

    pval.adjusted <- p.adjust(p.unadjusted, method = adjust.methods)
    pval <- min(pval.adjusted)
    max.index <- which(p.unadjusted == min(p.unadjusted), arr.ind = TRUE)

  }


  if(max & adjust.methods == "asymp"){

    #Calculating the covariace matrix

    tst.rslt1 <- rbind(tst.rslt[1,],subset(tst.rslt, grepl("FH", tst.rslt$W)))

    Z.tst.rslt1 <- tst.rslt1$Z
    q.tst.rslt1 <- tst.rslt1$Q
    var.tst.rslt1 <- tst.rslt1$Var

    wt1 <- c(list(a0=c(0,0)), wt)
    combo.wt <- combn(wt1,2)

    combo.wt.list <- list()
    for(i in 1:ncol(combo.wt)){combo.wt.list[[i]] <- combo.wt[,i]}

    combo.wt.list.up <- lapply(combo.wt.list,function(a){mapply('+',a)})



    wt2 <- lapply(combo.wt.list.up, function(a){apply(a,1,'sum')/2})
    d1 <- data.frame(do.call(rbind,wt2))


    wt3 <- unique(wt2)
    d2 <- data.frame(do.call(rbind,wt3))

    fit2<- ten(survfit(Surv(time, status) ~ arm, data = data.anal))

    #Testing (for calculating the covariances)

    # sink("out2.txt")
    comp(fit2, p= sapply(wt3, function(x){x[1]}), q= sapply(wt3, function(x){x[2]}))
    # sink()

    tst.rsltt <- attr(fit2, 'lrt')
    tst.rslt2 <- subset(tst.rsltt, grepl("FH", tst.rsltt$W))

    cov.tst.rslt11 <- tst.rslt2$Var
    d2$V <- cov.tst.rslt11


    d1d2 <- full_join(d1,d2, by = c("X1","X2"))

    cov.tst.rslt1 <- d1d2$V

    cov.tst.1 <- matrix(NA, nrow=length(wt1), ncol=length(wt1))


    cov.tst.1[lower.tri(cov.tst.1, diag=FALSE)] <- cov.tst.rslt1
    cov.tst <- t(cov.tst.1)
    cov.tst[lower.tri(cov.tst, diag=FALSE)] <- cov.tst.rslt1

    diag(cov.tst) <- var.tst.rslt1
    cov.tst.1 <- matrix(nearPD(cov.tst)$mat, length(Z.tst.rslt1),length(Z.tst.rslt1))
    #print(cov.tst.1)

    #z.val <- as.vector(ginv(Re(sqrtm(cov.tst)))%*%q.tst.rslt1)

    z.max <- max(abs(tst.rslt1$Z))
    cor.tst <- cov2cor(cov.tst.1)
    #print(cor.tst)

    #p.value=P(min(Z) < min(z.val))= 1 - P(Z_i >= min(z.val); for all i)

    if(one.sided){pval2 <- 1 - pmvnorm(lower = rep(-z.max, length(Z.tst.rslt1)),
                                       upper = rep(z.max, length(Z.tst.rslt1)),
                                       corr = cor.tst,
                                       algorithm = GenzBretz(maxpts=50000,abseps=0.00001))[1]
    max.tst <- which(abs(Z.tst.rslt1) == max(abs(Z.tst.rslt1)), arr.ind = TRUE)

    if(Z.tst.rslt1[max.tst] >= 0){pval <- 1 - pval2/2}
    if(Z.tst.rslt1[max.tst] < 0){pval <- pval2/2}

    }

    if(!one.sided){pval <- 1 - pmvnorm(lower = rep(-z.max, length(Z.tst.rslt1)),
                                       upper= rep(z.max, length(Z.tst.rslt1)),
                                       corr= cor.tst,
                                       algorithm= GenzBretz(maxpts=50000,abseps=0.00001))[1]}

    p.unadjusted <- stats::pnorm(q=tst.rslt1$Z)
    max.index <- which(p.unadjusted == min(p.unadjusted), arr.ind = TRUE)
    #max.index <- which(abs(Z.tst.rslt1) == max(abs(Z.tst.rslt1)), arr.ind = TRUE)

  }

  #Weighted log-rank test (FH weight): only one weight

  if(!max){

    tst.rslt1 <- subset(tst.rslt, grepl("FH", tst.rslt$W))
    pval <- stats::pnorm(q=tst.rslt1$Z)
    p.unadjusted <- pval
    max.index <- NULL

  }

  #Estimation (average HR)

  wt.rslt <- data.frame(attr(fit, 'lrw'))


  if(max){

    col.wt <- which(grepl("FH", colnames(wt.rslt))==TRUE, arr.ind = T)
    wt.rslt1.1 <- wt.rslt[,c(1,col.wt)]
    wt.rslt1 <- wt.rslt1.1[,max.index]

  }


  if(!max){
    col.wt <- which(grepl("FH", colnames(wt.rslt))==TRUE, arr.ind = T)
    wt.rslt1 <- wt.rslt[,col.wt]
  }

  #Performing weighted cox


  data.anal.event <- subset(data.anal, status==1)
  data.anal.cens <- subset(data.anal, status==0)

  data.anal.event <- data.anal.event[order(data.anal.event$time),]

  #Handing ties and matching length of weights

  wt.u <- unlist(wt.rslt1)


  time.freq <- as.matrix(table(data.anal.event$time))

  wt.un <- data.frame(cbind(as.numeric(row.names(time.freq)), time.freq[,1], wt.u))
  wt.all <- wt.un[rep(1:nrow(wt.un), times=wt.un[,2]),]
  data.anal.event$wt <- wt.all$wt.u

  #data.anal.event$wt <- unlist(wt.rslt1)


  data.anal.cens$wt <- rep(1, nrow(data.anal.cens))

  data.anal.w <- rbind(data.anal.event, data.anal.cens)
  data.anal.w$wt[data.anal.w$wt==0] <- 0.000001
  data.anal.w$wt2 <- -log(data.anal.w$wt)
  data.anal.w$id <- 1:nrow(data.anal.w)

  FH.est <- coxph(Surv(time, status) ~ arm + offset(wt2)+ cluster(id), weights=wt,
                  method= ties.method, data = data.anal.w)

  hr <- summary(FH.est)
  hr.est <- hr$conf.int[1]
  hr.est.se <- hr$coefficients[4]
  hr.low <- hr$conf.int[3]
  hr.up <- hr$conf.int[4]


  if(max){

    #Bonferroni adjustment

    hr.low.adjusted.BF <- exp(log(hr.est) - (stats::qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[4])
    hr.up.adjusted.BF <-  exp(log(hr.est) + (stats::qnorm(1- (alpha/(length(wt) + 1))))*hr$coefficients[4])

  }

  #Simultaneous CI using Hoetling T^2 and Exact MVN

  if(max & adjust.methods == "asymp"){

    set.seed(1234)

    c.star <- round(qmvnorm(1- alpha, corr=cor.tst, tail = "lower.tail")$quantile ,2)

    # EXACT   SCI

    hr.low.adjusted.E <- exp(log(hr.est) - c.star*hr.est.se)
    hr.up.adjusted.E  <- exp(log(hr.est) + c.star*hr.est.se)

  }

  hr.est.HT <- NULL
  hr.low.HT <- NULL
  hr.up.HT <- NULL

  if(HT.est & max & adjust.methods == "asymp"){

    nmodel <- length(wt) + 1
    prob.selection <- NULL
    hr.est2 <- NULL
    hr.est2.se <- NULL
    all.index <- 1:nmodel
    means <- rep(0,length(all.index))

    for(i in 1:nmodel){

      #Probability of model i is selected
      all.cindex <- all.index[-i]
      B <- matrix(rep(0,nmodel*nmodel),nrow=nmodel)
      diag(B) <- rep(1,nmodel)
      B[,i] <- -1
      B <- B[-i,]
      Bmeans <- as.vector(B %*% means)
      Bcov <-  B%*%cov.tst.1%*% t(B)
      Bcor <- cov2cor(Bcov)
      prob.selection[i] <- pmvnorm(lower=rep(0, length(all.cindex)),
                                   upper=Inf,
                                   corr=Bcor,
                                   algorithm = GenzBretz(maxpts=50000,abseps=0.00001))

      wt.rslt2 <- wt.rslt1.1[,i]

      data.anal.event2 <- subset(data.anal, status==1)
      data.anal.cens2 <- subset(data.anal, status==0)

      data.anal.event2 <- data.anal.event2[order(data.anal.event2$time),]

      #Handing ties and matching length of weights

      wt.u2 <- unlist(wt.rslt2)


      time.freq2 <- as.matrix(table(data.anal.event2$time))

      wt.un2 <- data.frame(cbind(as.numeric(row.names(time.freq2)), time.freq2[,1], wt.u2))
      wt.all2 <- wt.un2[rep(1:nrow(wt.un2), times=wt.un2[,2]),]
      data.anal.event2$wt <- wt.all2$wt.u2


      data.anal.cens2$wt <- rep(1, nrow(data.anal.cens2))

      data.anal.w2 <- rbind(data.anal.event2, data.anal.cens2)
      data.anal.w2$wt[data.anal.w2$wt==0] <- 0.000001
      data.anal.w2$wt2 <- -log(data.anal.w2$wt)


      FH.est2 <- coxph(Surv(time, status) ~ arm + offset(wt2), weights=wt,  method= ties.method, data = data.anal.w2)

      hr2 <- summary(FH.est2)
      hr.est2[i] <- hr2$conf.int[1]
      hr.est2.se[i] <- hr2$coefficients[3]


    }


    prob.selection2 <- prob.selection/sum(prob.selection)
    hr.est.HT <- exp(as.numeric(prob.selection2%*%log(hr.est2)))
    hr.est.HT.se <- sqrt(as.numeric(prob.selection2%*%hr.est2.se^2) + as.numeric(prob.selection2%*%log(hr.est2)^2) - (log(hr.est.HT))^2)
    hr.low.HT <- exp(log(hr.est.HT) - (stats::qnorm(1- alpha)*hr.est.HT.se))
    hr.up.HT <-  exp(log(hr.est.HT) + (stats::qnorm(1- alpha)*hr.est.HT.se))

  }

  if(max & adjust.methods != "asymp"){out <- list(pval=pval, pval.adjusted = pval.adjusted, p.unadjusted = p.unadjusted, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted.BF= hr.low.adjusted.BF, hr.up.adjusted.BF=hr.up.adjusted.BF,  Z.tst.rslt1= Z.tst.rslt1)}
  if(HT.est & max & adjust.methods == "asymp"){out <- list(cor=cor.tst,pval=pval, p.unadjusted= p.unadjusted, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted.BF= hr.low.adjusted.BF, hr.up.adjusted.BF=hr.up.adjusted.BF,  Z.tst.rslt1= Z.tst.rslt1, q.tst.rslt1=q.tst.rslt1, max.abs.z=z.max, hr.est.HT=hr.est.HT, hr.low.HT= hr.low.HT, hr.up.HT=hr.up.HT, prob.selection=prob.selection, data.anal.w=data.anal.w, hr.low.adjusted.E=hr.low.adjusted.E, hr.up.adjusted.E= hr.up.adjusted.E)}
  if(!HT.est & max & adjust.methods == "asymp"){out <- list(cor=cor.tst,pval=pval, p.unadjusted= p.unadjusted, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up, max.index=max.index, hr.low.adjusted.BF= hr.low.adjusted.BF, hr.up.adjusted.BF=hr.up.adjusted.BF,  Z.tst.rslt1= Z.tst.rslt1, q.tst.rslt1=q.tst.rslt1, max.abs.z=z.max, data.anal.w=data.anal.w, hr.low.adjusted.E=hr.low.adjusted.E, hr.up.adjusted.E= hr.up.adjusted.E)}
  if(!max){out <- list(pval=pval, hr.est=hr.est, hr.low=hr.low, hr.up=hr.up)}

  return(out)

}

sim_gsd_wlr <- function(N = ceiling(383.2),
                        enrollRates = tibble::tibble(Stratum = "All", duration = 12, rate = N/12),
                        failRates = tibble::tibble(Stratum = "All",
                                                   duration = c(4, 100),
                                                   failRate = log(2) / 15,  # median survival 15 month
                                                   hr = c(1, .6),
                                                   dropoutRate = 0.001),
                        rg = tibble(rho = 0, gamma = 0),
                        time = c(12, 24, 36),
                        lower = c(-0.694584174323271, 1.00239973891686, 1.99297019852853),
                        upper = c(3.7103028732625, 2.51140703831069, 1.99297019852853)
){


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

    # Cox model
    fit_cox <- coxph(Surv(time = tte, event = event)~ Treatment + strata(Stratum), data = sim_cut)
    cox_coef <- fit_cox$coefficients

    # Combo WLR
    res = rm.combo.WLRmax(time = sim_cut$tte,
                          status = sim_cut$event,
                          arm = sim_cut$Treatment,
                          wt = list(a1 = c(rg$rho, rg$gamma)),
                          max = FALSE)

    # Weighted log rank test
    z <- sim_cut %>% tensurv(txval = "Experimental") %>% tenFH(rg = rg)

    bind_cols(n = N, t = t, d = d, z, ahr = res$hr.est, ahr_cox = exp(cox_coef))
  }

  res <- bind_rows(lapply(time, foo, sim = sim))
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

path <- "/SFS/user/ctc/zhanyilo/gsdmvn"
load(file.path(path, "simulation/simu_gsd_wlr_boundary.Rdata") )

result <- list()
for(i in seq_along(x_weight)){
  x <- x_weight[[i]]
  N <- ceiling(max(x$N))
  lower <- x$Z[x$Bound == "Lower"]
  upper <- x$Z[x$Bound == "Upper"]
  wlr_fit <- sim_gsd_wlr(N = N, rg = tibble(rho = rho[i], gamma = gamma[i]),
                         lower = lower, upper = upper)
  result[[i]] <- wlr_fit

}

names(result) <- paste0("s", 1:length(result))

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(result, file = filename)


#----------------------
# HPC Submission code
#----------------------

# cd /SFS/scratch/zhanyilo/boundary
# rm *
# cp ~/gsdmvn/simulation/*.R .
# module add R/4.0.2
# qsub -t 1:2000 ~/runr.sh simu_gsd_wlr_boundary.R


#-------------------------------
# Summarize simulation results
#-------------------------------
#
# path <- "/SFS/scratch/zhanyilo/boundary"
#
# res <- list()
# for(i in 1:10000){
#   load(file.path(path, paste0(i, ".Rdata")))
#   try(
#     res[[i]] <- bind_rows(result, .id = "scenario")
#   )
# }
#
# res <- bind_rows(res)
# res %>% group_by(scenario, n, t, rho, gamma) %>%
#   summarise(events   = mean(d),
#             ahr      = mean(ahr),
#             lower    = mean(lower),
#             upper    = mean(upper)) %>%
#   group_by(scenario, n, rho, gamma) %>%
#   mutate(lower = cumsum(lower),
#          upper = cumsum(upper)) %>% data.frame() %>%
#   mutate_if(is.numeric, round, digits = 2)
#
#
# save(rho, gamma, x_weight, simu_res, file = "simulation/simu_gsd_wlr_boundary.Rdata")
