#' Fixed design sample size
#'
#' Computes fixed design sample size for many sample size methods.
#' Returns a `tibble` with a basic summary
#' @param x Sample size method; default is "AHR"; see examples and details.
#' @param alpha One-sided Type I error (strictly between 0 and 1)
#' @param power Power (`NULL` to compute power or strictly between 0 and `1 - alpha` otherwise)
#' @param ratio Experimental:Control randomization ratio
#' @param studyDuration study duration
#' @param ...
#' 
#' @return
#' @export
#'
#' @examples
#' # Average hazard ratio
#' y <- fixed_design("AHR", 
#'                   alpha = .025, Power = .9, 
#'                   enrollRates = tibble::tibble(Stratum = "All",  duration = 18, rate = 1),
#'                   failRates = tibble::tibble(Stratum = "All", duration = c(4, 100), failRate = log(2) / 12, hr = c(1, .6), dropoutRate = .001),
#'                   studyDuration = 36)
#' y %>% summary_fix()            
#' 
#' # Lachin and Foulkes (uses gsDesign::nSurv())
#' y <- fixed_design("LF", 
#'                   alpha = .025, Power = .9, 
#'                   enrollRates = tibble::tibble(Stratum = "All",  duration = 18, rate = 1),
#'                   failRates = tibble::tibble(Stratum = "All", duration = 100, failRate = log(2) / 12, hr = .7, dropoutRate = .001),
#'                   studyDuration = 36)
#' 
fixed_design <- function(x = c("AHR", "FH", "MB", "LF", "RD"), 
                         alpha = 0.025, power = NULL, ratio = 1, studyDuration = 36, ...){
   # --------------------------------------------- #
   #     check inputs                              #
   # --------------------------------------------- #
   if(!methods::hasArg(enrollRates)){
      stop("fixed_design: please input enrollRates!")
   }
   if(!methods::hasArg(failRates)){
      stop("fixed_design: please input failRates!")
   }
   
   x <- match.arg(x)
   
   y <- switch(x, 
               "AHR" = {
                  if (!is.null(power)){
                     d <- gs_design_ahr(alpha = alpha, beta = 1 - power,
                                        upar = qnorm(1 - alpha), lpar = -Inf,
                                        enrollRates = enrollRates,
                                        failRates = failRates,
                                        ratio  = ratio, 
                                        analysisTimes = studyDuration)
                  }else{
                     d <- gs_power_ahr(upar = qnorm(1 - alpha), lpar = -Inf,
                                       enrollRates = enrollRates,
                                       failRates = failRates,
                                       ratio  = ratio, 
                                       analysisTimes = studyDuration,
                                       events = NULL)
                  }
                  ans <- tibble::tibble(Design = "AHR",
                                        N = (d$analysis %>% filter(hypothesis == "H0"))$N,
                                        Events = (d$analysis %>% filter(hypothesis == "H0"))$Events,
                                        Time = (d$analysis %>% filter(hypothesis == "H0"))$Time,
                                        Bound = (d$bounds %>% filter(Bound == "Upper" & hypothesis == "H1"))$Z,
                                        alpha = (d$bounds %>% filter(hypothesis == "H0", Bound == "Upper"))$Probability,
                                        Power = (d$bounds %>% filter(hypothesis == "H1", Bound == "Upper"))$Probability)
                  
                  list(enrollRates = d$enrollRates, failRates = d$failRates, analysis = ans, design = "AHR")
                  },
               
               "FH" = {
                  temp1 <- methods::hasArg(weight)
                  temp2 <- methods::hasArg(rho)
                  temp3 <- methods::hasArg(gamma)
                  if(temp2 + temp3 == 0 & temp1 == 0){
                     weight <- function(x, arm0, arm1){gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)}
                  }
                  if(temp2 + temp3 >=1 & temp1 == 0){
                     weight <- function(x, arm0, arm1){gsdmvn:::wlr_weight_fh(x, arm0, arm1, 
                                                                              rho = ifelse(methods::hasArg(rho), rho, 0), 
                                                                              gamma =ifelse(methods::hasArg(gamma), gamma, 0.5))}
                  }
                  if (!is.null(power)){
                     d <- gs_design_wlr(alpha = alpha, beta = 1 - power,
                                        upar = qnorm(1 - alpha), lpar = -Inf,
                                        enrollRates = enrollRates, 
                                        failRates = failRates,
                                        ratio = ratio, 
                                        weight = weight,
                                        analysisTimes = studyDuration)
                  }else{
                     d <- gs_power_wlr(upar = qnorm(1 - alpha), lpar = -Inf,
                                       enrollRates = enrollRates, 
                                       failRates = failRates,
                                       ratio = ratio, 
                                       weight = weight,
                                       analysisTimes = studyDuration,
                                       events = NULL)
                  }
                  ans <- tibble::tibble(Design = "FH",
                                        N = (d$analysis %>% filter(hypothesis == "H0"))$N,
                                        Events = (d$analysis %>% filter(hypothesis == "H0"))$Events,
                                        Time = (d$analysis %>% filter(hypothesis == "H0"))$Time,
                                        Bound = (d$bounds %>% filter(Bound == "Upper" & hypothesis == "H1"))$Z,
                                        alpha = (d$bounds %>% filter(hypothesis == "H0", Bound == "Upper"))$Probability,
                                        Power = (d$bounds %>% filter(hypothesis == "H1", Bound == "Upper"))$Probability)
                  
                  list(enrollRates = d$enrollRates, failRates = d$failRates, analysis = ans, design = "FH")
                  },
               
               
               "MB" = {
                  # check if power is NULL or not
                  if(!is.null(power)){
                     d <- gs_design_wlr(alpha = alpha,
                                        beta = 1 - power,
                                        enrollRates = enrollRates, 
                                        failRates = failRates,
                                        ratio = 1, 
                                        weight = function(x, arm0, arm1){
                                           gsdmvn:::wlr_weight_fh(x, arm0, arm1, 
                                                                  rho = ifelse(methods::hasArg(rho), rho, 0),
                                                                  gamma = ifelse(methods::hasArg(gamma), gamma, 0),
                                                                  tau = ifelse(methods::hasArg(tau), tau, 6))},
                                        upper = gs_b,
                                        upar = qnorm(1 - alpha),
                                        lower = gs_b,
                                        lpar = -Inf,
                                        analysisTimes = studyDuration) 
                  }else{
                     d <- gs_power_wlr(enrollRates = enrollRates, 
                                       failRates = failRates,
                                       ratio = 1, 
                                       weight = function(x, arm0, arm1){
                                          gsdmvn:::wlr_weight_fh(x, arm0, arm1, 
                                                                 rho = ifelse(methods::hasArg(rho), rho, 0),
                                                                 gamma = ifelse(methods::hasArg(gamma), gamma, 0),
                                                                 tau = ifelse(methods::hasArg(tau), tau, 6))},
                                       upper = gs_b,
                                       upar = qnorm(1 - alpha),
                                       lower = gs_b,
                                       lpar = -Inf,
                                       analysisTimes = studyDuration,
                                       events = NULL) 
                  }
                  # get the output of MB
                  ans <- tibble::tibble(Design = "MB",
                                        N = (d$analysis %>% filter(hypothesis == "H0"))$N,
                                        Events = (d$analysis %>% filter(hypothesis == "H0"))$Events,
                                        Time = (d$analysis %>% filter(hypothesis == "H0"))$Time,
                                        Bound = (d$bounds %>% filter(Bound == "Upper" & hypothesis == "H1"))$Z,
                                        alpha = (d$bounds %>% filter(hypothesis == "H0", Bound == "Upper"))$Probability,
                                        Power = (d$bounds %>% filter(hypothesis == "H1", Bound == "Upper"))$Probability)
                  
                  list(enrollRates = d$enrollRates, failRates = d$failRates, analysis = ans, design = "MB")
                  
                  
               },
               
                 
               "LF" = { # Should do checks of inputs here
                  # calculate the S: durations of piecewise constant event rates 
                  m <- length(failRates$failRate)
                  if (m == 1){S <- NULL}else{S <- failRates$duration[1:(m-1)]}
                  
                  # calculate the ahr as the hr in nSurv
                  if (!is.null(power)){
                     dd <- gs_design_ahr(alpha = alpha, beta = 1 - power,
                                         upar = qnorm(1 - alpha), lpar = -Inf,
                                         enrollRates = enrollRates,
                                         failRates = failRates,
                                         ratio  = ratio, 
                                         analysisTimes = studyDuration)
                  }else{
                     dd <- gs_power_ahr(upar = qnorm(1 - alpha), lpar = -Inf,
                                        enrollRates = enrollRates,
                                        failRates = failRates,
                                        ratio  = ratio, 
                                        analysisTimes = studyDuration,
                                        events = NULL)
                  }
                  
                  # use nSuve to develop the design
                  d <- gsDesign::nSurv(alpha = alpha, beta = if(is.null(power)){NULL}else{1 - power}, 
                                       ratio = ratio, hr = dd$analysis$AHR[1],
                                       # failRates
                                       lambdaC = failRates$failRate,
                                       S = S, eta = failRates$dropoutRate,  
                                       # enrollRates
                                       gamma = enrollRates$rate, R = enrollRates$duration,
                                       T = studyDuration, minfup = studyDuration - sum(enrollRates$duration))
                  
                  ans <- tibble::tibble(Design = "LF",
                                        N = d$n,
                                        Events = d$d,
                                        Time = d$T,
                                        alpha = d$alpha,
                                        Power = d$power)
                  list(enrollRates = d$enrollRates, failRates = d$failRates, analysis = ans, design = "LF")
                  },
               
               
               "RD" = {
                  list(sum_ = tibble::tibble(Option = 5, Design = "RD") # everything that needs to be returned should be in the tibble (p_C, p_E, N, alpha, and power)
                )},
                 
               
               list(sum_= tibble::tibble(Option = 99, Design = "Not implemented"),
                                           x = "Enter implemented design type in x"))
   class(y) <- "fixed_design"
   return(y)    
}

# summary function of the fixed design
summary <- function(y){
   return(y$analysis)
}