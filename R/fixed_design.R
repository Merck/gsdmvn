#' Fixed design sample size
#'
#' Computes fixed design sample size for many sample size methods.
#' Returns a `tibble` with a basic summary
#' @param x Sample size method; default is "AHR"; see examples and details.
#' @param alpha One-sided Type I error (strictly between 0 and 1)
#' @param power Power (`NULL` to compute power or strictly between 0 and `1 - alpha` otherwise)
#' @param ratio Experimental:Control randomization ratio
#' @param studyDuration study duration
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
fixed_design <- function(x = "AHR", alpha = 0.025, power = .9, ratio = 1, studyDuration = 36, ...){
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
                    
                     list(sum_ = (d$analysis %>% filter(hypothesis == "H0"))$N,
                          enrollRates = d$enrollRates)
                     },
                 
                 "FH" = { # This will call gs_design_wlr or gs_power_wlr
                    list(sum_ = tibble::tibble(Option = 2, Design = "FH"),
                         c = 2)
                    },
                 
                 "MB" = {
                    list(sum_ = tibble::tibble(Option = 3, Design = "MB"),
                         stuff = "stuff")
                    },
                 
                 
                 "LF" = { # Should do checks of inputs here
                    m <- length(failRates$failRate)
                    if (m == 1){S <- NULL}else{S <- failRates$duration[1:(m-1)]}
                    
                    d <- gsDesign::nSurv(
                       alpha = alpha, 
                       beta = ifelse(is.null(power), NULL, 1 - power), 
                       hr = failRates$hr,  
                       S = S, 
                       lambdaC = failRates$failRate, 
                       eta = failRates$dropoutRate,
                       R = enrollRates$duration, 
                       gamma = enrollRates$rate,
                       T = studyDuration, 
                       ratio = ratio)
                    
                    list(sum_ = tibble::tibble(Option = 4, Design = "LF"),
                         x = d)
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