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
NULL
#' Grid points for group sequential design numerical integration
#'
#' Points and weights for Simpson's rule numerical integration from
#' p 349 - 350 of Jennison and Turnbull book.
#' This is not used for arbitrary integration, but for the canonical form of Jennison and Turnbull.
#' mu is computed elsewhere as drift parameter times sqrt of information.
#' Since this is a lower-level routine, no checking of input is done; calling routines should
#' ensure that input is correct.
#' Lower limit of integration can be \code{-Inf} and upper limit of integration can be \code{Inf}
#'
#' @details
#' Jennison and Turnbull (p 350) claim accuracy of \code{10E-6} with \code{r=16}.
#' The numerical integration grid spreads out at the tail to enable accurate tail probability calcuations.
#'
#'
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param mu Mean of normal distribution (scalar) under consideration
#' @param a lower limit of integration (scalar)
#' @param b upper limit of integration (scalar \code{> a})
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Define odd numbered grid points for real line.
#'    \item Trim points outside of [a, b] and include those points.
#'    \item If extreme, include only 1 point where density will be essentially 0.
#'    \item Define even numbered grid points between the odd ones.
#'    \item Compute weights for odd numbered grid points.
#'    \item Combine odd- and even-numbered grid points with their corresponding weights.
#'    \item Return a tibble of with grid points in z and numerical integration weights in z.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @return A \code{tibble} with grid points in \code{z} and numerical integration weights in \code{w}
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # approximate variance of standard normal (i.e., 1)
#' gridpts_() %>% summarise(var = sum(z^2 * w * dnorm(z)))
#'
#' # approximate probability above .95 quantile (i.e., .05)
#' gridpts_(a = qnorm(.95), b = Inf) %>% summarise(p05 = sum(w * dnorm(z)))
gridpts_ <- function(r = 18, mu = 0, a = -Inf, b = Inf){
  # Define odd numbered grid points for real line
  x <- c(mu - 3 - 4 * log(r / (1:(r - 1))),
         mu - 3 + 3 * (0:(4 * r)) / 2 / r,
         mu + 3 + 4 * log(r / (r - 1):1)
         )
  # Trim points outside of [a, b] and include those points
  if (min(x) < a) x <- c(a, x[x > a])
  if (max(x) > b) x <- c(x[x < b], b)
  # If extreme, include only 1 point where density will be essentially 0
  m <- length(x)
  if (m == 1) return(tibble::tibble(z=x, w=1))

  # Define even numbered grid points between the odd ones
  y <- (x[2:m] + x[1:(m-1)]) / 2

  # Compute weights for odd numbered grid points
  i <- 2:(m-1)
  wodd <- c(x[2] - x[1],
            (x[i + 1] - x[i - 1]),
            x[m] - x[m - 1]) / 6

  weven <- 4 * (x[2:m] - x[1:(m-1)]) / 6

  # Now combine odd- and even-numbered grid points with their
  # corresponding weights
  z <- rep(0, 2*m - 1)
  z[2 * (1:m) - 1] <- x
  z[2 * (1:(m-1))] <- y
  w <- z
  w[2 * (1:m) - 1] <- wodd
  w[2 * (1:(m-1))] <- weven

  return(tibble::tibble(z=z, w=w))
}
