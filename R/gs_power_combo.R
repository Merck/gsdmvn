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

#' Group sequential design power using MaxCombo test under non-proportional hazards
#'
#' @inheritParams gs_design_combo
#' @inheritParams pmvnorm_combo
#'
#' @examples
#' library(dplyr)
#' library(mvtnorm)
#'
#' enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 500/12)
#'
#' failRates <- tibble::tibble(Stratum = "All",
#'                             duration = c(4, 100),
#'                             failRate = log(2) / 15,  # median survival 15 month
#'                             hr = c(1, .6),
#'                             dropoutRate = 0.001)
#'
#' fh_test <- rbind( data.frame(rho = 0, gamma = 0, tau = -1,
#'                              test = 1,
#'                              Analysis = 1:3,
#'                              analysisTimes = c(12, 24, 36)),
#'                   data.frame(rho = c(0, 0.5), gamma = 0.5, tau = -1,
#'                              test = 2:3,
#'                              Analysis = 3, analysisTimes = 36)
#' )
#'
#' # User defined bound
#' gs_power_combo(enrollRates, failRates, fh_test, upar = c(3,2,1), lpar = c(-1, 0, 1))
#'
#' # Minimal Information Fraction derived bound
#' gs_power_combo(enrollRates, failRates, fh_test,
#'                upper = gs_spending_combo,
#'                upar  = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
#'                lower = gs_spending_combo,
#'                lpar  = list(sf = gsDesign::sfLDOF, total_spend = 0.2))
#'
#' @importFrom mvtnorm GenzBretz
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate if lower and upper bounds have been specified.
#'    \item Extract info, info_fh, theta_fh and corr_fh from utility.
#'    \item Extract sample size via the maximum sample size of info.
#'    \item Calculate information fraction either for fixed or group sequential design.
#'    \item Compute spending function using \code{gs_bound()}.
#'    \item Compute probability of crossing bounds under the null and alternative
#'     hypotheses using \code{gs_prob_combo()}.
#'    \item Export required information for boundary and crossing probability
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @export
gs_power_combo <- function(enrollRates,
                           failRates,
                           fh_test,
                           ratio = 1,
                           binding = FALSE,
                           upper = gs_b,
                           upar = c(3,2,1),
                           lower = gs_b,
                           lpar = c(-1, 0, 1),
                           algorithm = GenzBretz(maxpts= 1e5, abseps= 1e-5),
                           ...){

  # Currently only support user defined lower and upper bound
  stopifnot( identical(upper, gs_b) | identical(upper, gs_spending_combo) )
  stopifnot( identical(lower, gs_b) | identical(lower, gs_spending_combo) )

  # Obtain utilities
  utility <- gs_utility_combo(enrollRates = enrollRates,
                              failRates = failRates,
                              fh_test = fh_test,
                              ratio = ratio,
                              algorithm = algorithm, ...)

  info     <- utility$info_all
  info_fh  <- utility$info
  theta_fh <- utility$theta
  corr_fh  <- utility$corr

  # Sample size
  n <- max(info$N)

  # Information Fraction
  if(length(unique(fh_test$Analysis)) == 1){
    # Fixed design
    min_info_frac <- 1

  }else{

    info_frac <- tapply(info$info0, info$test, function(x) x / max(x))
    min_info_frac <- apply(do.call(rbind, info_frac), 2, min)
  }

  # Obtain spending function
  bound <- gs_bound(alpha = upper(upar, min_info_frac),
                    beta = lower(lpar, min_info_frac),
                    analysis = info_fh$Analysis,
                    theta = theta_fh * sqrt(n),
                    corr = corr_fh,
                    binding_lower_bound = binding,
                    algorithm = algorithm,
                    alpha_bound = identical(upper, gs_b),
                    beta_bound = identical(lower, gs_b),
                    ...)


  # Probability Cross Boundary uner Alternative
  prob <- gs_prob_combo(upper_bound = bound$upper,
                        lower_bound = bound$lower,
                        analysis = info_fh$Analysis,
                        theta = theta_fh * sqrt(n),
                        corr = corr_fh,
                        algorithm = algorithm, ...)

  # Probability Cross Boundary under Null
  prob_null <- gs_prob_combo(upper_bound = bound$upper,
                             lower_bound = if(binding){bound$lower}else{rep(-Inf, nrow(bound))},
                             analysis = info_fh$Analysis,
                             theta = rep(0, nrow(info_fh)),
                             corr = corr_fh,
                             algorithm = algorithm, ...)

  if(binding == FALSE){
    prob_null$Probability[prob_null$Bound == "Lower"] <- NA
  }

  prob$Probability_Null <- prob_null$Probability

  # Prepare output
  db <- merge(data.frame(Analysis = 1:(nrow(prob)/2), prob, Z = unlist(bound)),
              unique(info_fh[, c("Analysis", "Time", "N", "Events")])
  )

  db[order(db$Bound, decreasing = TRUE), c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability", "Probability_Null")]

}
