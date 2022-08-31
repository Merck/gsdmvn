#  Copyright (c) 2022 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.
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
#' @importFrom stats qnorm uniroot
#' @importFrom utils tail
#' @importFrom gsDesign gsDesign sfLDOF
#' @importFrom gsDesign2 eAccrual AHR tEvents
#' @importFrom dplyr select left_join n

NULL
#' Fixed and group sequential design under non-proportional hazards
#'
#' \code{gs_design_nph()} is a flexible routine to provide a sample size or power for a fixed or
#' group sequential design under
#' various non-proportional hazards assumptions for either single or multiple strata studies.
#' The piecewise exponential distribution allows a simple method to specify a distribtuion
#' and enrollment pattern
#' where the enrollment, failure and dropout rates changes over time.
#' @param enrollRates Piecewise constant enrollment rates by stratum and time period.
#' @param failRates Piecewise constant control group failure rates, duration for each piecewise constant period,
#' hazard ratio for experimental vs control, and dropout rates by stratum and time period.
#' @param ratio Experimental:Control randomization ratio
#' @param alpha One-sided Type I error
#' @param beta Type II error
#' @param analysisTimes Final calendar time if beta is not NULL
#' @param IF information fraction at planned analyses
#' @param upper Function to produce upper bound
#' @param lower Function to produce lower bound
#' @param upar Parameters to pass to upper
#' @param lpar Parameters to pass to lower
#' @param r  Control for grid size; normally leave at default of \code{r=18}
#'
#' @return A list with 3 tibbles: 1) \code{enrollRates} with \code{enrollRates$rate} adjusted by sample size calculation and adding \code{N} with
#' cumulative enrollment at end of each enrollment rate period,
#' 2) \code{failRates} as input,
#' 3) code{bounds} with a row for each bound and each analysis; rows contain the variables \code{Analysis} with analysis number, \code{Z} with Z-value bound,
#' \code{Probability} with cumulative probability of crossing bound at each analysis under the alternate hypothesis input,
#' \code{theta} standardized effect size at each analysis, \code{info} cumulative statistical information for \code{theta} at each analysis,
#' \code{Time} expected timing of analysis,
#' \code{avehr} expected average hazard ratio at time of analysis,
#' \code{Events} expected events an time of analysis under alternate hypothesis,
#' \code{info0} information under null hypothesis with same expected total events under alternate hypothesis, and
#' \code{N} expected enrollment at time of analysis.
#' @examples
#' library(dplyr)
#' 
#' # Design
#' library(dplyr)
#' 
#' x <- gs_design_nph()
#' # Notice cumulative power is 0.9 (90%)
#' x
#' # Check Type I error, non-binding; should be .025 (2.5%)
#' gs_power_nph(enrollRates = x$enrollRates,
#'            failRates = x$failRates %>% mutate(hr = 1),
#'            events = (x$bounds %>% filter(Bound == "Upper"))$Events,
#'            upar = (x$bounds %>% filter(Bound == "Upper"))$Z,
#'            lpar = rep(-Inf,3),
#'            upper = gs_b,
#'            lower = gs_b
#'            ) %>% filter(abs(Z) < Inf)
#' # Power under proportional hazards, HR = 0.75
#' gs_power_nph(enrollRates = x$enrollRates,
#'            failRates = x$failRates %>% mutate(hr = .75),
#'            events = (x$bounds %>% filter(Bound == "Upper"))$Events,
#'            upar = (x$bounds %>% filter(Bound == "Upper"))$Z,
#'            lpar = rep(-Inf,3),
#'            upper = gs_b,
#'            lower = gs_b
#'            ) %>% filter(abs(Z) < Inf)
#'
#' @export
#'
gs_design_nph <- function(enrollRates=tibble::tibble(Stratum="All",
                                                   duration=c(2,2,10),
                                                   rate=c(3,6,9)),
                        failRates=tibble::tibble(Stratum="All",
                                                 duration=c(3,100),
                                                 failRate=log(2)/c(9,18),
                                                 hr=c(.9,.6),
                                                 dropoutRate=rep(.001,2)),
                        ratio = 1, # NOT YET IMPLEMENTED
                        alpha = 0.025,       # One-sided Type I error
                        beta = 0.1,          # NULL if enrollment is not adapted
                        analysisTimes = 30,  # CURRENTLY GIVES ONLY FINAL CALENDAR TIME OF TRIAL
                        IF = c(.25, .75, 1), # relative information fraction timing (vector, if not NULL)
                        upper = gs_b,
                        # Default is Lan-DeMets approximation of
                        upar = gsDesign::gsDesign(k=3, test.type=1, n.I=c(.25, .75, 1), maxn.IPlan = 1,
                                                  sfu=sfLDOF, sfupar = NULL)$upper$bound,
                        lower = gs_b,
                        lpar = c(qnorm(.1), rep(-Inf, length(IF) - 1)),
                        r = 18
){
  errbeta <- function(x = 1, info, theta, Zupper, Zlower, upar, lpar, beta){
     1 -  beta - max((gs_prob(theta = theta, upper = Zupper, lower = Zlower, upar = upar, lpar = lpar,
                             info = x * info, r = r) %>% filter(Bound == "Upper"))$Probability)
  }
  avehr <- gsDesign2::AHR(enrollRates = enrollRates,
                          failRates = failRates,
                          totalDuration = analysisTimes,
                          ratio = ratio)
  K <- max(length(IF), length(analysisTimes))
  finalEvents <- max(avehr$Events)
  if (is.null(IF)){
    avehr$IF <- avehr$Events/finalEvents
  }else{
    avehr <- NULL
    for(i in seq_along(IF)){
      avehr <- rbind(avehr,
                     gsDesign2::tEvents(enrollRates = enrollRates,
                                        failRates = failRates,
                                        ratio = ratio,
                                        targetEvents = IF[i] * finalEvents)
               )
    }
    avehr$IF <- IF
  }
  targ <- (qnorm(alpha) + qnorm(beta))^2 / log(avehr$AHR[nrow(avehr)])^2 * (1 + ratio)^2 / ratio
  interval <- c(.9, 1.5) * targ / finalEvents

# Now we can solve for the inflation factor for the enrollment rate to achieve the desired power
  res <- try(
    uniroot(errbeta,
            interval = interval,
            theta = -log(avehr$AHR),
            Zupper = upper,
            Zlower = lower,
            upar = upar,
            lpar = lpar,
            info = avehr$info,
            beta = beta,
            tol = .0001
    )
  )
  if(inherits(res,"try-error")){stop("gs_design_nph(): Sample size solution not found")}
  enrollRates = enrollRates %>% mutate(rate = rate * res$root)
  avehr <- gsDesign2::AHR(enrollRates,
                          failRates = failRates,
                          totalDuration = avehr$Time) %>%
           mutate(Analysis = 1:n())
  N <- cumsum(enrollRates$rate * enrollRates$duration) %>% tail(1)
  # Once eAccrual is fixed, should not need pmin in the following
  avehr$N <- pmin(N, gsDesign2::eAccrual(avehr$Time, enrollRates))
  bounds <-
    gs_prob(theta = -log(avehr$AHR),
            upper = upper,
            lower = lower,
            upar = upar,
            lpar = lpar,
            info = avehr$info
    )
  return(list(enrollRates=enrollRates %>% mutate(N=cumsum(duration * rate)),
              failRates=failRates,
              bounds = bounds %>% left_join(avehr %>% select(-info), by="Analysis"))
  )
}
