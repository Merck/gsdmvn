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

#' MaxCombo Group sequential boundary crossing probabilities
#'
#' @inheritParams pmvnorm_combo
#' @param upper_bound a numeric vector of upper bound
#' @param lower_bound a numeric vector of lower bound
#' @param analysis an integer vector of the interim analysis index
#' @param theta a numeric vector of effect size under alternative hypothesis
#' @param corr a matrix of correlation matrix
#'
#' @importFrom mvtnorm GenzBretz
#'
gs_prob_combo <- function(upper_bound,
                          lower_bound,
                          analysis,
                          theta,
                          corr,
                          algorithm = GenzBretz(maxpts= 1e5, abseps= 1e-5),
                          ...){

  n_analysis <- length(unique(analysis))

  p <- c()
  q <- c()
  for(k in 1:n_analysis){
    k_ind <- analysis <= k


    # Efficacy Bound
    if(k == 1){
      lower <- upper_bound[1]
      upper <- Inf
    }else{
      lower <- c(lower_bound[1:(k-1)], upper_bound[k])
      upper <- c(upper_bound[1:(k-1)], Inf)
    }


    p[k] <- pmvnorm_combo(lower,
                          upper,
                          group = analysis[k_ind],
                          mean = theta[k_ind],
                          corr = corr[k_ind, k_ind])

    # Futility Bound
    if(k == 1){
      lower <- -Inf
      upper <- lower_bound[k]
    }else{
      lower <- c(lower_bound[1:(k-1)], -Inf)
      upper <- c(upper_bound[1:(k-1)], lower_bound[k])
    }

    q[k] <- pmvnorm_combo(lower,
                          upper,
                          group = analysis[k_ind],
                          mean  = theta[k_ind],
                          corr  = corr[k_ind, k_ind])

  }

  data.frame(Bound = rep(c("Upper", "Lower"), each = n_analysis),
             Probability = c(cumsum(p),cumsum(q)))

}


