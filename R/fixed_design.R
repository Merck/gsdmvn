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
fixed_design <- function(x = "AHR", alpha = 0.025, power = .9, ratio = NULL, studyDuration = NULL, ...){
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
                                       analysisTimes = studyDuration)
                  }
                  ans <- tibble::tibble(Design = "FH",
                                        N = (d$analysis %>% filter(hypothesis == "H0"))$N %>% round(1),
                                        Events = (d$analysis %>% filter(hypothesis == "H0"))$Events %>% round(1),
                                        Time = (d$analysis %>% filter(hypothesis == "H0"))$Time %>% round(2),
                                        Bound = (d$bounds %>% filter(Bound == "Upper" & hypothesis == "H1"))$Z %>% round(4),
                                        alpha = (d$bounds %>% filter(hypothesis == "H0"))$Probability %>% round(4),
                                        Power = (d$bounds %>% filter(hypothesis == "H1"))$Probability %>% round(4))
                  
                  list(analysis = ans, enrollRates = d$enrollRates, failRates = d$failRates, design = "AHR")
                  },
               
               "FH" = { # This will call gs_design_wlr or gs_power_wlr
                  # if the weight is not inputted, set it as the default
                  if(!methods::hasArg(weight)){
                     weight = function(x, arm0, arm1){gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)}
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
                                       analysisTimes = studyDuration)
                  }
                  ans <- tibble::tibble(Design = "FH",
                                        N = (d$analysis %>% filter(hypothesis == "H0"))$N %>% round(1),
                                        Events = (d$analysis %>% filter(hypothesis == "H0"))$Events %>% round(1),
                                        Time = (d$analysis %>% filter(hypothesis == "H0"))$Time %>% round(2),
                                        Bound = (d$bounds %>% filter(Bound == "Upper" & hypothesis == "H1"))$Z %>% round(4),
                                        alpha = (d$bounds %>% filter(hypothesis == "H0"))$Probability %>% round(4),
                                        Power = (d$bounds %>% filter(hypothesis == "H1"))$Probability %>% round(4))
                  
                  list(analysis = ans, enrollRates = d$enrollRates, failRates = d$failRates, design = "FH")
                  },
               
               
               "MB" = {
                  list(sum_ = tibble::tibble(Option = 3, Design = "MB"), stuff = "stuff")
               },
               
                 
               "LF" = { # Should do checks of inputs here
                  m <- length(failRates$failRate)
                  if (m == 1){S <- NULL}else{S <- failRates$duration[1:(m-1)]}
                  
                  d <- gsDesign::nSurv(alpha = alpha, beta = ifelse(is.null(power), NULL, 1 - beta), 
                                       ratio = ratio,
                                       # failRates
                                       hr = failRates$hr,  lambdaC = failRates$failRate,
                                       S = S, eta = failRates$dropoutRate,  
                                       # enrollRates
                                       gamma = enrollRates$rate, R = enrollRates$duration,
                                       T = studyDuration, minfup = studyDuration - sum(enrollRates$duration))
                  
                  ans <- tibble::tibble(Design = "LF",
                                        N = d$n %>% round(1),
                                        Events = d$d %>% round(1),
                                        Time = d$T %>% round(2),
                                        alpha = d$alpha %>% round(4),
                                        Power = d$power %>% round(4))
                  list(analysis = ans, enrollRates = d$enrollRates, failRates = d$failRates, design = "FH")
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
# let start with the non-S3 method first. Since we need to align how can S3 class we need later.
summary_fix <- function(y){
   y$sum_
}