#' Sample Size Calculation based on RMST method 
#' 
#' @examples 
#' # set enrollment rates
#' enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 500/12)
#' 
#' # set failure rates
#' failRates <- tibble::tibble(
#'   Stratum = "All",
#'   duration = c(4, 100),
#'   failRate = log(2) / 15,  # median survival 15 month
#'   hr = c(1, .6),
#'   dropoutRate = 0.001)
#'
#' fixed_design_size_rmst(enrollRates, failRates, analysisTimes = 36)
#' fixed_design_power_rmst(enrollRates, failRates, analysisTimes = 36)
#' fixed_design_size_rmst(enrollRates, failRates, analysisTimes = 36, beta = 1 - 0.887)
fixed_design_size_rmst <- function(enrollRates, 
                                   failRates, 
                                   analysisTimes,
                                   ratio = 1,
                                   alpha = 0.025, 
                                   beta = 0.1){
  
  gs_arm <- gs_create_arm(enrollRates, failRates, 
                          ratio = ratio,              # Randomization ratio
                          total_time = analysisTimes) # Total study duration
  
  arm0 <- gs_arm[["arm0"]]
  arm1 <- gs_arm[["arm1"]]
  
  # Sample size for RMST at cut point 36. 
  npsurv <- npsurvSS::size_two_arm(arm0, arm1, 
                                   power = 1 - beta, 
                                   alpha = alpha, 
                                   test = list(test="rmst difference", 
                                               milestone = arm0$total_time)) 
  
  analysis <- tibble::tibble(
    design = "RMST", 
    N = npsurv[["n"]], 
    Events = npsurv[["d"]], 
    Time = analysisTimes, 
    Bound = - qnorm(alpha),
    alpha = alpha, 
    Power = 1 - beta
  )
  
  res <- list(enrollRates = enrollRates, 
              failRates = failRates, 
              analysis = analysis, 
              design = "RMST")
  class(res) <- c("fixed_design", "list")
  
  res
}



fixed_design_power_rmst <- function(enrollRates, 
                                    failRates, 
                                    analysisTimes,
                                    ratio = 1,
                                    alpha = 0.025){
  
  gs_arm <- gs_create_arm(enrollRates, failRates, 
                          ratio = ratio,              # Randomization ratio
                          total_time = analysisTimes) # Total study duration
  
  n <- sum(enrollRates$duration * enrollRates$rate)
  n0 <- n / (ratio + 1)
  n1 <- n - n0

  
  arm0 <- gs_arm[["arm0"]]
  arm1 <- gs_arm[["arm1"]]
  arm0$size <- n0 
  arm1$size <- n1
  
  d <- prob_event.arm(arm0, tmax =  arm0$total_time) * n0 + 
    prob_event.arm(arm1, tmax =  arm0$total_time) * n1
  
  # Sample size for RMST at cut point 36. 
  npsurv <- npsurvSS::power_two_arm(arm0, arm1, 
                                    alpha = alpha, 
                                    test = list(test="rmst difference", 
                                                milestone = arm0$total_time)) 
  
  analysis <- tibble::tibble(
    design = "RMST", 
    N = n, 
    Events = d, 
    Time = analysisTimes, 
    Bound = - qnorm(alpha),
    alpha = alpha, 
    Power = npsurv
  )
  
  res <- list(enrollRates = enrollRates, 
              failRates = failRates, 
              analysis = analysis, 
              design = "RMST")
  class(res) <- c("fixed_design", "list")
  
  
  res
}