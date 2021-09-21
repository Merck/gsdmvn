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

#' @importFrom stats pnorm
#' @importFrom dplyr "%>%" summarise filter
#' @importFrom tibble tibble
NULL
#' Group sequential boundary crossing probabilities
#'
#' @param theta natural parameter for group sequentia design representing expected drift at time of each analysis
#' @param upper function to compute upper bound
#' @param lower function to compare lower bound
#' @param upar parameter to pass to upper
#' @param lpar parameter to pass to lower
#' @param info statistical information at each analysis
#' @param r  Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @return A `tibble` with a row for each finite bound and analysis containing the following variables:
#' Analysis analysis number
#' Bound Upper (efficacy) or Lower (futility)
#' Z Z-value at bound
#' Probability probability that this is the first bound crossed under the given input
#' theta approximate natural parameter value required to cross the bound
#'
#' @details
#' Approximation for \code{theta} is based on Wald test and assumes the observed information is equal to the expected.
#' @examples
#' # Asymmetric 2-sided design
#' gs_prob(theta = 0, upar = rep(2.2, 3), lpar = rep(0, 3), 
#'         upper=gs_b, lower=gs_b,  info = 1:3)
#' # One-sided design
#' x <- gs_prob(theta = 0, upar = rep(2.2, 3), lpar = rep(-Inf, 3), 
#'              upper=gs_b, lower=gs_b,  info = 1:3)
#' # Without filtering, this shows unneeded lower bound
#' x
#' # Filter to just show bounds intended for use
#' x %>% filter(abs(Z) < Inf)
#' @export
gs_prob <- function(theta, upper=gs_b, lower=gs_b, upar, lpar, info, r = 18){
   # deal with R cmd check messages
   Z <- h <- NULL
   K <- length(info)
   Zupper <- upper(upar, info)
   Zlower <- lower(lpar, info)
   if (length(theta) == 1) theta <- rep(theta, K)
   upperProb <- rep(NA, K)
   lowerProb <- rep(NA, K)
   for(k in seq_along(info)){
     if(k==1){
        upperProb[1] <- if(Zupper[1] < Inf) {pnorm(Zupper[1], mean = sqrt(info[1]) * theta[1], lower.tail = FALSE)}else{0}
        lowerProb[1] <- if(Zlower[1] > -Inf){pnorm(Zlower[1], mean = sqrt(info[1]) * theta[1])}else{0}
        g <- h1(r = r, theta = theta[1], I = info[1], a = Zlower[1], b = Zupper[1])
     }else{
        # Cross upper bound
        upperProb[k] <- if(Zupper[k]< Inf){
              hupdate(r = r, theta = theta[k], I = info[k], a = Zupper[k], b = Inf,
                             thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = g) %>%
              summarise(sum(h)) %>% as.numeric()
           }else{0}
      # Cross lower bound
        lowerProb[k] <- if(Zlower[k] > -Inf){
              hupdate(r = r, theta = theta[k], I = info[k], a = -Inf, b = Zlower[k],
                             thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = g) %>%
              summarise(sum(h)) %>% as.numeric()
        }else{0}
        # if k < K, update numerical integration for next analy
        if (k < K) g <- hupdate(r = r, theta = theta[k], I = info[k], a = Zlower[k], b = Zupper[k],
                                       thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = g)
     }
   }
     return(tibble::tibble(
              Analysis = rep(1:K, 2),
              Bound = c(rep("Upper", K), rep("Lower", K)),
              Z= c(Zupper, Zlower),
              Probability = c(cumsum(upperProb),
                              cumsum(lowerProb)),
              theta = rep(theta, 2),
              info = rep(info, 2)) # %>% filter(abs(Z) < Inf)
         )
}
