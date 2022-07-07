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

#' @importFrom stats dnorm
#' @importFrom tibble tibble
NULL
#' Update numerical integration for group sequential design
#'
#' Update grid points for numerical integration from one analysis to the next
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param theta Drift parameter for current analysis
#' @param I Information at current analysis
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#' @param thetam1  Drift parameter for previous analysis
#' @param Im1 Information at previous analysis
#' @param gm1 numerical integration grid from \code{h1()} or previous run of \code{hupdate()}
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Compute the square root of the change in information.
#'    \item Compute the grid points for group sequential design numerical integration.
#'    \item Update the integration.
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
#' # 2nd analysis with no interim bound and drift 0 should have mean 0, variance 1
#' hupdate_() %>% summarise(mu = sum(z * h), var = sum((z - mu)^2 * h))
hupdate_ <- function(r = 18, theta = 0, I = 2, a = -Inf, b = Inf, thetam1 = 0, Im1 = 1, gm1 = h1_()){
  # sqrt of change in information
  rtdelta <- sqrt(I - Im1)
  rtI <- sqrt(I)
  rtIm1 <- sqrt(Im1)
  g <- gridpts_(r = r, mu = theta * rtI, a= a, b = b)
  # update integration
  mu <- theta * I - thetam1 * Im1
  h <- rep(0, length(g$z))
  for(i in seq_along(g$z)){
    x <- (g$z[i] * rtI - gm1$z * rtIm1 - mu ) / rtdelta
    h[i] <- sum(gm1$h * dnorm(x))
  }
  h <- h * g$w * rtI / rtdelta
  return(tibble::tibble(z = g$z, w = g$w, h = h))
}
