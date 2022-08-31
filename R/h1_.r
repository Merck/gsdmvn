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

#' @importFrom stats dnorm pnorm
#' @importFrom tibble tibble
NULL
#' Initialize numerical integration for group sequential design
#'
#' Compute grid points for first interim analysis in a group sequential design
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param theta Drift parameter for first analysis
#' @param I Information at first analysis
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#'
#' @details Mean for standard normal distribution under consideration is \code{mu = theta * sqrt(I)}
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Compute drift at analysis 1.
#'    \item Compute deviation from drift.
#'    \item Compute standard normal density, multiply by grid weight.
#'    \item Return a tibble of z, w, and h.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @return A \code{tibble} with grid points in \code{z}, numerical integration weights in \code{w},
#' and a normal density with mean \code{mu = theta * sqrt{I}} and variance 1 times the weight in \code{w}.
#'
#' @examples
#' library(dplyr)
#' # Replicate variance of 1, mean of 35
#' gsdmvn:::h1_(theta = 5, I = 49) %>% summarise(mu = sum(z * h), var = sum((z - mu)^2 * h))
#'
#' # Replicate p-value of .0001 by numerical integration of tail
#' gsdmvn:::h1_(a = qnorm(.9999)) %>% summarise(p = sum(h))
h1_ <- function(r = 18, theta = 0, I = 1, a = -Inf, b = Inf){
  # fix for binding message
  z <- w <- h <- NULL
  # compute drift at analysis 1
  mu <- theta * sqrt(I);
  g <- gridpts_(r, mu, a, b)
  # compute deviation from drift
  x <- g$z - mu
  # compute standard normal density, multiply by grid weight and return
  # values needed for numerical integration
  return(tibble::tibble(z = g$z, w = g$w, h = g$w * dnorm(x)))
}
