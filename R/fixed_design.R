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
fixed_design <- function(x = "AHR", alpha = 0.025, power = NULL, ratio = NULL, studyDuration = NULL, ...){
   # --------------------------------------------- #
   #     check inputs                              #
   # --------------------------------------------- #
   if(!methods::hasArg(enrollRates)){
      stop("fixed_design: please input enrollRates!")
   }
   if(!methods::hasArg(failRates)){
      stop("fixed_design: please input failRates!")
   }
   
   if(is.null(ratio)){
      ratio <- 1
      message("fixed_design: the default 1:1 randomization ratio is used!")
   }
   if(is.null(studyDuration)){
      studyDuration <- 36
      message("fixed_design: the default 36 months study duration is used!")
   }
   
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
               
               "FH" = { # This will call gs_design_wlr or gs_power_wlr
                  # if the weight is not inputted, set it as the default
                  if(!methods::hasArg(weight)){
                     weight <- function(x, arm0, arm1){gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)}
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
                  # check if we have rho, gamma, and tau input
                  if(!methods::hasArg(rho)){
                     message("fixed_design: rho is not input and the default value rho = 0 is used!")
                  }
                  if(!methods::hasArg(gamma)){
                     message("fixed_design: gamma is not input and the default value gamma = 0 is used!")
                  }
                  if(!methods::hasArg(tau)){
                     message("fixed_design: tau is not input and the default value tau = 6 is used!")
                  }
                  
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
summary_fix <- function(y){
   return(y$analysis)
}