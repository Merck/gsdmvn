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

#' @importFrom tibble tibble
#' @importFrom dplyr lag
NULL
#' Information and effect size based on risk difference for binary outcome
#'
#' Based on binomial success or failure rates and sample size for control and experimental arms,
#' statistical information under the null and alternate hypothesis are computed. 
#' Supports stratified population with information weighted by either statistical information (default) or 
#' population size. 
#' Multiple sample sizes can be input for a group sequential design.
#' Method matches Fleiss, Tytun and Ury for superiority and Farrington and Manning for 
#' non-inferiority or super-superiority evaluations of a single stratum.
#' Extends these methods by allowing stratified population and alternate weighting approaches to strata.
#' @param p_control probability of an event in the control group under alternate hypothesis; between 0 and 1; 
#' can be a vector with different strata. 
#' @param p_experimental probability of an event in the experimental group under alternate hypothesis; between 0 and 1; 
#' can be a vector with different strata.
#' @param theta0 real value between xx and xx; null hypothesis for `p_control - p_experimental`; natural parameter with default = 0.
#' While one or both of `p_control` and `theta0` can vary between strata; can be a single value or one per stratum. 
#' @param ratio Experimental:Control randomization ratio; real value > 0.
#' @param N vector of increasing positive real values with total sample size at each analysis; not required to be integer since formulas are asymptotic.
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate that input `N` is a real vector with increasing values > 0.
#'    \item Validate separately that `p_control` and `p_experimental` have length 1 or the same length as `N`; each value should be between 0 and 1.
#'    either constant or changing over time. If changing, this should reflect average event rate cumulatively
#'    through each analysis.
#'    \item Validate that `ratio` is a single real value > 0.
#'    \item If `p_control` or `p_experimental` is length 1, make each a vector with same length as `N` by repeating values.
#'    \item Compute `theta = p_control - p_experimental`, the natural parameter under the alternate hypothesis.
#'    \item Compute proportion of observations assigned to control and save in `xi`; 
#'    this is assumed to be common across strata.
#'    \item At each analysis, compute null and alternate hypothesis rates for control and experimental arms by stratum:
#'    \itemize{
#'       \item Compute 
#'    }
#'    \item Return a tibble of Analysis, N, "Risk difference", theta, info, info0.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @return a \code{tibble} with columns \code{Analysis, N, "Risk difference", p_control, p_experimental,
#' theta, info, info0.}
#' \code{info, info0} contains statistical information under H1, H0, respectively.
#' "Risk difference" For analysis \code{k}, \code{Time[k]} is the maximum of \code{analysisTimes[k]} and the expected time
#' required to accrue the targeted \code{events[k]}.
#' \code{AHR} is expected average hazard ratio at each analysis.
#' @details The \code{AHR()} function computes statistical information at targeted event times.
#' The \code{tEvents()} function is used to get events and average HR at targeted \code{analysisTimes}.
#' @export
#'
#' @examples
#' library(gsDesign)
#' library(gsDesign2)
#' 
#' # Only put in targeted events
#' gs_info_ahr(events = c(30, 40, 50))
#' # Only put in targeted analysis times
#' gs_info_ahr(analysisTimes = c(18, 27, 36))
#' # Some analysis times after time at which targeted events accrue
#' # Check that both Time >= input analysisTime and Events >= input events
#' gs_info_ahr(events = c(30, 40, 50), analysisTimes = c(16, 19, 26))
#' gs_info_ahr(events = c(30, 40, 50), analysisTimes = c(14, 20, 24))
gs_info_ahr <- function(enrollRates=tibble::tibble(Stratum="All",
                                                  duration=c(2,2,10),
                                                  rate=c(3,6,9)),
                       failRates=tibble::tibble(Stratum="All",
                                                duration=c(3,100),
                                                failRate=log(2)/c(9,18),
                                                hr=c(.9,.6),
                                                dropoutRate=rep(.001,2)),
                       ratio=1,               # Experimental:Control randomization ratio
                       events = NULL,         # Events at analyses
                       analysisTimes = NULL   # Times of analyses
){
  ################################################################################
  # Check input values
  K <- 0
  if (is.null(analysisTimes) && is.null(events)) stop("gs_info_ahr(): One of
                                                      events and analysisTimes must be a
                                                      numeric value or vector with increasing values")
  if (!is.null(analysisTimes)){
    if (!is.numeric(analysisTimes) || !is.vector(analysisTimes) || min(analysisTimes - dplyr::lag(analysisTimes, def=0))<=0
       )stop("gs_info_ahr(): analysisTimes must be NULL a numeric vector with positive increasing values")
    K <- length(analysisTimes)
  }
  if (!is.null(events)){
    if (!is.numeric(events) || !is.vector(events) || min(events - dplyr::lag(events, default=0))<=0
    )stop("gs_info_ahr(): events must be NULL or a numeric vector with positive increasing values")
    if(K==0){
      K <- length(events)
    }else if (K != length(events)) stop("gs_info_ahr(): If both events and analysisTimes specified, must have same length")
  }
  # end check input values
  ################################################################################
  avehr <- NULL
  if(!is.null(analysisTimes)){
    avehr <- gsDesign2::AHR(enrollRates = enrollRates, failRates = failRates, ratio = ratio,
                            totalDuration = analysisTimes)
    for(i in seq_along(events)){
      if (avehr$Events[i] < events[i]){
        avehr[i,] <- gsDesign2::tEvents(enrollRates = enrollRates, failRates = failRates, ratio = ratio,
                                 targetEvents = events[i])
      }
    }
  }else{
     for(i in seq_along(events)){
        avehr <- rbind(avehr,
                   gsDesign2::tEvents(enrollRates = enrollRates, failRates = failRates, ratio = ratio,
                                      targetEvents = events[i]))
     }
  }
  avehr$Analysis <- 1:nrow(avehr)
  avehr$theta = -log(avehr$AHR)
  return(avehr %>% dplyr::transmute(Analysis, Time, Events, AHR, theta, info, info0))
}
