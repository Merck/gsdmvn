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

#' Create "npsurvSS" arm object
#'
#' @param total_time total analysis time
#' @inheritParams gs_info_ahr
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate if there is only one stratum.
#'    \item Calculate the accrual duration.
#'    \item calculate the accrual intervals.
#'    \item Calculate the accrual parameter as the proportion of enrollment rate*duration.
#'    \item Set cure proportion to zero.
#'    \item set survival intervals and shape.
#'    \item Set fail rate in failRates to the Weibull scale parameter for the survival distribution in the arm 0.
#'    \item Set the multiplication of hazard ratio and fail rate to the Weibull scale parameter
#'    for the survival distribution in the arm 1.
#'    \item Set the shape parameter to one as the exponential distribution for
#'    shape parameter for the loss to follow-up distribution
#'    \item Set the scale parameter to one as the scale parameter for the loss to follow-up
#'     distribution since the exponential distribution is supported only
#'    \item Create arm 0 using \code{npsurvSS::create_arm()} using the parameters for arm 0.
#'    \item Create arm 1 using \code{npsurvSS::create_arm()} using the parameters for arm 1.
#'    \item Set the class of the two arms.
#'    \item Return a list of the two arms.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
gs_create_arm <- function(enrollRates,
                          failRates,
                          ratio,
                          total_time = 1e6){

  n_stratum <- length(unique(enrollRates$Stratum))
  if(n_stratum > 1){
    stop("Only one stratum is supported")
  }

  accr_time     <- sum(enrollRates$duration)
  accr_interval <- cumsum(enrollRates$duration)
  accr_param    <- enrollRates$rate * enrollRates$duration / sum(enrollRates$rate * enrollRates$duration)

  surv_cure     <- 0                    # No cure proportion
  surv_interval <- c(0, c(utils::head(failRates$duration, -1), Inf))
  surv_shape    <- 1                    # Exponential Distribution
  surv_scale0   <- failRates$failRate
  surv_scale1   <- failRates$hr * failRates$failRate

  loss_shape    <- 1                         # Exponential Distribution
  loss_scale    <- failRates$dropoutRate[1]  # Only Exponential Distribution is supported

  # Control Group
  arm0 <- npsurvSS::create_arm(size = 1,

                               accr_time = accr_time,
                               accr_dist = "pieceuni",
                               accr_interval = accr_interval,
                               accr_param = accr_param,

                               surv_cure = surv_cure,
                               surv_interval = surv_interval,
                               surv_shape = surv_shape,
                               surv_scale = surv_scale0,

                               loss_shape = loss_shape,
                               loss_scale = loss_scale,

                               total_time = total_time)


  # Control Group
  arm1 <- npsurvSS::create_arm(size = ratio,

                               accr_time = accr_time,
                               accr_dist = "pieceuni",
                               accr_interval = accr_interval,
                               accr_param = accr_param,

                               surv_cure = surv_cure,
                               surv_interval = surv_interval,
                               surv_shape = surv_shape,
                               surv_scale = surv_scale1,

                               loss_shape = loss_shape,
                               loss_scale = loss_scale,

                               total_time = total_time)

  class(arm0) <- c("list", "arm")
  class(arm1) <- c("list", "arm")

  list(arm0 = arm0,
       arm1 = arm1)

}

