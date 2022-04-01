#  Copyright (c) 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
#
#  This file is part of the gsdmvn program.
#
#  gsdmvn is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Group sequential design using weighted log-rank test under non-proportional hazards
#'
#' @inheritParams gs_design_ahr
#' @inheritParams gs_info_wlr
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate if input analysisTimes is a positive number or a positive increasing sequence.
#'    \item Validate if input IF is a positive number or positive increasing sequence on (0, 1] with final value of 1.
#'    \item Validate if inputs IF and analysisTimes  have the same length if both have length > 1.
#'    \item Compute information at input analysisTimes using \code{gs_info_wlr()}.
#'    \item Compute sample size and bounds using \code{gs_design_npe()}.
#'    \item Return a list of design enrollment, failure rates, and bounds.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @export
#' 
#' @examples
#' library(dplyr)
#' library(mvtnorm)
#' library(gsDesign)
#' 
#' enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 500/12)
#' 
#' failRates <- tibble::tibble(Stratum = "All",
#'                             duration = c(4, 100),
#'                             failRate = log(2) / 15,  # median survival 15 month
#'                             hr = c(1, .6),
#'                             dropoutRate = 0.001)
#' 
#' x <- gsDesign::gsSurv(
#'   k = 3, 
#'   test.type = 4, 
#'   alpha = 0.025, beta = 0.2, 
#'   astar = 0, timing = 1,
#'   sfu = sfLDOF, sfupar = 0, 
#'   sfl = sfLDOF, sflpar = 0, 
#'   lambdaC = 0.1, 
#'   hr = 0.6, hr0 = 1, 
#'   eta = 0.01, gamma = 10,
#'   R = 12, S = NULL,
#'   T = 36, minfup = 24, 
#'   ratio = 1)
#' 
#' # User defined boundary
#' gs_design_wlr(
#'   enrollRates = enrollRates, 
#'   failRates = failRates,
#'   ratio = 1, 
#'   alpha = 0.025, beta = 0.2,
#'   weight = function(x, arm0, arm1){
#'     gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)},
#'   upar = x$upper$bound,
#'   lpar = x$lower$bound,
#'   analysisTimes = c(12, 24, 36))
#' 
#' # Boundary derived by spending function
#' gs_design_wlr(
#'   enrollRates = enrollRates, 
#'   failRates = failRates,
#'   ratio = 1, 
#'   alpha = 0.025, beta = 0.2,
#'   weight = function(x, arm0, arm1){
#'     gsdmvn:::wlr_weight_fh(x, arm0, arm1, rho = 0, gamma = 0.5)
#'   },
#'   upper = gs_spending_bound,
#'   upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
#'   lower = gs_spending_bound,
#'   lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.2),
#'   analysisTimes = c(12, 24, 36))

gs_design_wlr <- function(
  enrollRates = tibble::tibble(Stratum = "All",
                               duration = c(2,2,10),
                               rate = c(3,6,9)),
  failRates = tibble::tibble(Stratum = "All",
                             duration = c(3,100),
                             failRate = log(2)/c(9,18),
                             hr = c(.9,.6),
                             dropoutRate = rep(.001,2)),
  ratio = 1,             # Experimental:Control randomization ratio
  weight = wlr_weight_fh,
  approx = "asymptotic",
  alpha = 0.025,         # One-sided Type I error
  beta = 0.1,            # NULL if enrollment is not adapted
  IF = NULL,             # relative information fraction timing (vector, if not NULL; increasing to 1)
  analysisTimes = 36,    # Targeted times of analysis or just planned study duration
  binding = FALSE,
  upper = gs_b,
  # Default is Lan-DeMets approximation of
  upar = gsDesign(k = 3, test.type = 1, n.I = c(.25, .75, 1), sfu = sfLDOF, sfupar = NULL)$upper$bound,
  lower = gs_b,
  lpar = c(qnorm(.1), -Inf, -Inf), # Futility only at IA1
  h1_spending = TRUE,
  test_upper = TRUE,
  test_lower = TRUE,
  r = 18,
  tol = 1e-6
){
  # --------------------------------------------- #
  #     check input values                        #
  # --------------------------------------------- #
  msg <- "gs_design_wlr(): analysisTimes must be a positive number or positive increasing sequence"
  if (!is.vector(analysisTimes,mode = "numeric")) stop(msg)
  if (min(analysisTimes - dplyr::lag(analysisTimes, def=0))<=0) stop(msg)
  msg <- "gs_design_wlr(): IF must be a positive number or positive increasing sequence on (0, 1] with final value of 1"
  if (is.null(IF)){IF <- 1}
  if (!is.vector(IF,mode = "numeric")) stop(msg)
  if (min(IF - dplyr::lag(IF, def=0))<=0) stop(msg)
  if (max(IF) != 1) stop(msg)
  msg <- "gs_design_wlr(): IF and analysisTimes must have the same length if both have length > 1"
  if ((length(analysisTimes)>1) & (length(IF) > 1) & (length(IF) != length(analysisTimes))) stop(msg)
  
  # --------------------------------------------- #
  #     get information at input analysisTimes    #
  # --------------------------------------------- #
  y <- gs_info_wlr(
    enrollRates, failRates, 
    ratio = ratio, events = NULL, 
    analysisTimes = analysisTimes,
    weight = weight, approx = approx)
  
  finalEvents <- y$Events[nrow(y)]
  IFalt <- y$Events / finalEvents
  
  # Check if IF needed for (any) IA timing
  K <- max(length(analysisTimes), length(IF))
  nextTime <- max(analysisTimes)
  
  if(length(IF) == 1){IF <- IFalt}else{
    IFindx <- IF[1:(K-1)]
    for(i in seq_along(IFindx)){
      if(length(IFalt) == 1){y <-
        rbind(
          gsDesign2::tEvents(
            enrollRates, failRates, 
            targetEvents = IF[K - i] * finalEvents, ratio = ratio,
            interval = c(.01, nextTime)) %>% mutate(theta=-log(AHR), Analysis=K-i),
          y)
      }else if (IF[K-i] > IFalt[K-i]) y[K - i,] <-
          gsDesign2::tEvents(
            enrollRates, failRates, 
            targetEvents = IF[K - i] * finalEvents, ratio = ratio,
            interval = c(.01, nextTime)) %>%
          dplyr::transmute(Analysis = K - i, Time, Events, AHR, theta=-log(AHR), info, info0)
      nextTime <- y$Time[K - i]
    }
  }
  
  y$Analysis <- 1:K
  y$N <- gsDesign2::eAccrual(x = y$Time, enrollRates = enrollRates)
  if(h1_spending){
    theta1 <- y$theta
    info1 <- y$info
  }else{
    theta1 <- 0
    info1 <- y$info0
  }
  
  # --------------------------------------------- #
  #     combine all the calculations              #
  # --------------------------------------------- #
  suppressMessages(
  allout <- gs_design_npe(
    theta = y$theta, theta1 = theta1,
    info = y$info, info0 = y$info0, info1 = info1,
    alpha = alpha, beta = beta, binding = binding,
    upper = upper, upar = upar, test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol) %>%
    
    # Add Time, Events, AHR, N from gs_info_ahr call above
    full_join(y %>% select(-c(info, info0, theta)), by = "Analysis") %>%
    
    select(c("Analysis", "Bound", "Time",
             "N", "Events", 
             "Z", "Probability",
             "AHR", "theta", 
             "info", "info0", "info1", "IF", "hypothesis")) %>%  
    arrange(desc(hypothesis), desc(Bound), Analysis)  
  )
  
  allout$Events <- allout$Events * allout$info[K] / y$info[K]
  allout$N <- allout$N * allout$info[K] / y$info[K]
  
  # add `~HR at bound`, `HR generic` and `Nominal p`
  allout <- allout %>% mutate(
    "~HR at bound" = gsDesign::zn2hr(z = Z, n = Events, ratio = ratio),
    "Nominal p" = pnorm(-Z)
  ) 
  # --------------------------------------------- #
  #     get bounds to output                      #
  # --------------------------------------------- #
  bounds <- allout %>% 
    select(all_of(c("Analysis", "Bound", "Probability", "hypothesis", "Z",
                    "~HR at bound", "Nominal p" )))
  # --------------------------------------------- #
  #     get analysis summary to output            #
  # --------------------------------------------- #
  analysis <- allout %>% 
    select(Analysis, Time, N, Events, AHR, theta, info, IF, hypothesis) %>% 
    unique()
  
  # --------------------------------------------- #
  #     return the output                         #
  # --------------------------------------------- #
  output <- list(
    enrollRates = enrollRates %>% mutate(rate = rate * allout$info[K] / y$info[K]),
    failRates = failRates,
    bounds = bounds,
    analysis = analysis)
  class(output) <- c("wlr", "gs_design", class(output))
  return(output)
  
}