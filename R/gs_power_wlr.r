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

#' Group sequential design power using weighted log rank test under non-proportional hazards
#'
#' @inheritParams gs_design_wlr
#' @inheritParams gs_power_ahr
#' @section Specification:
#' 
#' \if{latex}{
#'  \itemize{
#'    \item Compute information and effect size for Weighted Log-rank test using \code{gs_info_wlr()}.
#'    \item Compute group sequential bound computation with non-constant effect using \code{gs_power_npe()}.
#'    \item Combine information and effect size and power and return a
#'    tibble  with columns Analysis, Bound, Time, Events, Z, Probability, AHR,  theta, info, and info0.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @export
#' @examples 
#' traget_events <- c(30, 40, 50)
#' target_analysisTime <- c(10, 24, 30)
#' 
#' # calculate the power for targeted number of events
#' gs_power_wlr(
#'   events = traget_events, 
#'   analysisTimes = NULL,
#'   upper = gs_b,
#'   upar = gsDesign::gsDesign(k = length(traget_events), test.type = 1, n.I = traget_events, maxn.IPlan = max(traget_events), sfu = sfLDOF, sfupar = NULL)$upper$bound,
#'   lower = gs_b,           
#'   lpar = c(qnorm(.1), rep(-Inf, length(events) - 1))) 
#'   
#' # calculate the power for targeted analysis time
#' gs_power_wlr(
#'   events = NULL, 
#'   analysisTimes = target_analysisTime,
#'   upper = gs_b,
#'   upar = gsDesign::gsDesign(k = length(traget_events), test.type = 1, n.I = traget_events, maxn.IPlan = max(traget_events), sfu = sfLDOF, sfupar = NULL)$upper$bound,
#'   lower = gs_b,           
#'   lpar = c(qnorm(.1), rep(-Inf, length(events) - 1))) 
#' 
#' # calculate the power for targeted analysis time & number of events
#' gs_power_wlr(
#'   events = traget_events, 
#'   analysisTimes = target_analysisTime,
#'   upper = gs_b,
#'   upar = gsDesign::gsDesign(k = length(traget_events), test.type = 1, n.I = traget_events, maxn.IPlan = max(traget_events), sfu = sfLDOF, sfupar = NULL)$upper$bound,
#'   lower = gs_b,           
#'   lpar = c(qnorm(.1), rep(-Inf, length(events) - 1))) 
#'   
gs_power_wlr <- function(
  enrollRates = tibble::tibble(
    Stratum = "All",
    duration = c(2, 2, 10),
    rate = c(3, 6, 9)),
  failRates = tibble::tibble(
    Stratum = "All",
    duration = c(3,100),
    failRate = log(2)/c(9, 18),
    hr = c(.9,.6),
    dropoutRate = rep(.001, 2)),
  ratio = 1,               # Experimental:Control randomization ratio
  weight = wlr_weight_fh,
  approx = "asymptotic",
  events = c(30, 40, 50), # Targeted events of analysis
  analysisTimes = NULL,   # Targeted times of analysis
  binding = FALSE,
  upper = gs_b,           # Default is Lan-DeMets approximation of
  upar = gsDesign::gsDesign(k = 3, test.type = 1, n.I = c(30, 40, 50), maxn.IPlan = 50, sfu = sfLDOF, sfupar = NULL)$upper$bound,
  test_upper = TRUE,
  lower = gs_b,           # Futility only at IA1
  lpar = c(qnorm(.1), rep(-Inf, length(events) - 1)), 
  test_lower = TRUE,
  r = 18,
  tol = 1e-6
){
  
  # get the number of analysis
  K <- max(length(events), length(analysisTimes), na.rm = TRUE)
  
  # ---------------------------------------- #
  #    calculate the asymptotic variance     #
  #       and statistical information        #
  # ---------------------------------------- #
  x <- gs_info_wlr(
    enrollRates = enrollRates,
    failRates = failRates,
    ratio = ratio,
    events = events,
    weight = weight,
    analysisTimes = analysisTimes
  )
  
  # ---------------------------------------- #
  #  given the above statistical information #
  #         calculate the power              #
  # ---------------------------------------- #
  y <- gs_power_npe(
    theta = x$theta, 
    info = x$info, 
    info0 = x$info0, 
    binding = binding,
    upper = upper, 
    lower = lower, 
    upar = upar, 
    lpar= lpar,
    test_upper = test_upper, 
    test_lower = test_lower,
    r = r, 
    tol = tol)
  
  # ---------------------------------------- #
  #         organize the outputs             #
  # ---------------------------------------- #
  # summarize the bounds
  suppressMessages(
  bounds <- y %>%
    left_join(x %>% select(Analysis, Events)) %>% 
    mutate(
      `~HR at bound` = gsDesign::zn2hr(z = Z, n = Events, ratio = ratio),
      `Nominal p` = pnorm(-Z)
    ) %>% 
    select(Analysis, Bound, Probability, Z, `~HR at bound`, `Nominal p`)
  )
  
  # summarize the analysis
  suppressMessages(
  analysis <- x %>% 
    select(Analysis, Time, Events, AHR) %>% 
    mutate(N = gsDesign2::eAccrual(x = x$Time, enrollRates = enrollRates)) %>% 
    left_join(y) %>% 
    select(Analysis, Time, N, Events, AHR, theta, info, IF)
  )
  
  output <- list(
    enrollRates = enrollRates, 
    failRates = failRates,
    bounds = bounds,
    analysis = analysis)
  
  class(output) <- c("wlr", "gs_design", class(output))
  
  return(output)
}
