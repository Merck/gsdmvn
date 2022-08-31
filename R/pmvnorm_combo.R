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

#' Multivariate Normal Distribution for Multivariate Maximum Statistics
#'
#' Computes the distribution function of the multivariate normal distribution
#' with maximum statistics for arbitrary limits and correlation matrices
#'
#' @inheritParams mvtnorm::pmvnorm
#' @param group the vector of test statistics group.
#' @param ... additional parameters transfer to `mvtnorm::pmvnorm`
#' @details
#'
#' Let $Z = {Z_ij}$ be a multivariate normal distribution.
#' Here i is a group indicator and j is a within group statistics indicator.
#' Let G_i = max({Z_ij}) for all test within one group.
#' This program are calculating the probability
#'
#'   $$Pr( lower < max(G) < upper )$$
#'
#' @importFrom mvtnorm GenzBretz
#'
#' @export
pmvnorm_combo <- function(lower,
                          upper,
                          group,
                          mean,
                          corr,
                          algorithm = GenzBretz(maxpts= 1e5, abseps= 1e-5),
                          ...){

  # Number of test in each group
  n_test <- as.numeric(table(group))


  # Ensure positive definitive of the correlation matrix
  if(! corpcor::is.positive.definite(corr)){
    corr <- corpcor::make.positive.definite(corr)
    corr <- stats::cov2cor(corr)
  }

  # One dimension test
  if(length(mean) == 1){
    p <- pnorm(mean, lower) - pnorm(mean, upper)
    return(p)
  }

  # One test for all group or lower bound is -Inf.
  if(all(n_test == 1) | all(lower == -Inf) ){
    p <- mvtnorm::pmvnorm(lower = rep(lower, n_test),
                 upper = rep(upper, n_test),
                 mean = mean,
                 corr = corr,
                 sigma = NULL,
                 algorithm = algorithm,
                 ...)

    return(p)

  # General Algorithm
  }else{

    # Re-arrange test based on the order for number of test
    group <- as.numeric(factor(group, order(n_test)))

    mean <- mean[order(group)]
    corr <- corr[order(group), order(group)]
    group <- group[order(group)]

    n_test <- as.numeric(table(group))



    # Split by number of test == 1
    lower1 <- lower[n_test == 1]
    upper1 <- upper[n_test == 1]

    lower2 <- lower[n_test > 1]
    upper2 <- upper[n_test > 1]

    # Combination of all possible test
    k <- length(lower2)
    test_ind <- split(matrix(c(1,-1), nrow = k, ncol = 2, byrow = TRUE), 1:k)
    test_ind <- expand.grid(test_ind)
    test <- split(test_ind, 1:nrow(test_ind))

    p <- sapply(test, function(x){
      lower_bound <- rep(c(lower1, rep(-Inf, k)), n_test)
      upper_bound <- rep(c(upper1, ifelse(x == 1, upper2, lower2)), n_test)

      p_bound <- mvtnorm::pmvnorm(lower = lower_bound,
                         upper = upper_bound,
                         mean = mean,
                         corr = corr,
                         sigma = NULL,
                         algorithm = algorithm,
                         ...)

      prod(x) * p_bound

    })

    return(sum(p))

  }

}
