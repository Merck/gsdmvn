#' Fixed design sample size
#'
#' Computes fixed design sample size for many sample size methods.
#' Returns a `tibble` with a basic summary
#' @param x Sample size method; default is "AHR"; see examples and details.
#' @param alpha One-sided Type I error (strictly between 0 and 1)
#' @param Power Power (`NULL` to compute power or strictly between 0 and `1 - alpha` otherwise)
#' @param ratio Experimental:Control randomization ratio
#' @param ... Additional parameters needed; see examples.
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
#'                   analysisTimes = 36)
#' y %>% summary() %>% gt::gt()            
#' 
fixed_design <- 
  function(x = "AHR", alpha = 0.025, Power = .9, ratio = 1, ...){
     y <- switch(x, 
                 "AHR" = {
                      if (!is.null(Power)){
                         d <- gs_design_ahr(alpha = alpha, beta = 1 - Power,
                                            upar = qnorm(1 - alpha), lpar = -Inf,
                                            enrollRates = enrollRates,
                                            failRates = failRates,
                                            ratio  = ratio, 
                                            analysisTimes = analysisTimes,...)
                       }else{
                           d <- gs_power_ahr(upar = qnorm(1 - alpha), lpar = -Inf,
                                             enrollRates = enrollRates,
                                             failRates = failRates,
                                             ratio  = ratio, 
                                             analysisTimes = analysisTimes,...)
                         }
                     list(sum_ = d$analysis %>% filter(Bound == "Upper", hypothesis == "H0"),
                         enrollRates = d$enrollRates)
              },
            "FH" = {
                list(sum_ = tibble::tibble(Option = 2, Design = "FH"),
                     c = 2)
                },
            "MB" = {
                list(sum_ = tibble::tibble(Option = 3, Design = "MB"),
                     stuff = "stuff")
              },
            "LF" = { 
                list(sum_ = tibble::tibble(Option = 4, Design = "LF"), 
                     magic = "LF")
              },
            "RD" = {
                list(sum_ = tibble::tibble(Option = 5, Design = "RD") # everything that needs to be returned should be in the tibble (p_C, p_E, N, alpha, and power)
                    )},
                list(sum_= tibble::tibble(Option = 99, Design = "Not implemented"),
                                          x = "Enter implemented design type in x")
                    )
     class(y) <- "fixed_design"
     return(y)    
   }
summary.fixed_design <- function(y){y$sum_}