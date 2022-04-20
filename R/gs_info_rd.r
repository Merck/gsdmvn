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
#' Information and effect size under risk difference 
#' 
#' @param p_c rate at the control group
#' @param p_e rate at the experimental group 
#' @param n sample size 
#' @param theta0 the standardized treatment effect under H0 
#' @param delta0 the treatment effect under H0 
#' @param ratio Experimental:Control randomization ratio
#'
#' @export
#' @examples
#' gs_info_rd()
#' 
gs_info_rd <- function(
  p_c = .15,
  p_e = .13,
  N = c(50, 80, 100),
  theta0 = 0,
  delta0 = 0, 
  ratio = 1  
){
  # set the sample size of control/experiment group
  n_e <- N / (1 + ratio)
  n_c <- N * ratio / (1 + ratio)
  # set the pooled rate
  p_pool <- p_c * n_c / N + p_e * n_e / N
  # set d
  d <- ifelse(p_c > p_e, 1, 1)
  # set the treatment effect 
  theta <- d * rep(p_c - p_e, length(N))
  
  # if it is superiority, non-inferiority design
  if(theta0 == 0){
    sigma2_H1 <- p_c * (1 - p_c) / n_c  + p_e * (1 - p_e) / n_e
    sigma2_H0 <- p_pool * (1 - p_pool) * (1 / n_c + 1 / n_e)
    info1 <- 1 / sigma2_H1
    info0 <- 1 / sigma2_H0
  }
  
  # if it is super-superiority designs
  if(theta0 > 0){
    p_e0 <- (p_c + ratio * p_e - d * theta0) / (ratio + 1)
    p_c0 <- p_e0 + d * theta0
    sigma2_H0 <- p_c0 * (1 - p_c0) / n_c  + p_e0 * (1 - p_e0) / n_e
    sigma2_H1 <- p_c * (1 - p_c) / n_c  + p_e * (1 - p_e) / n_e
    info1 <- 1 / sigma2_H1
    info0 <- 1 / sigma2_H0
  }
  
  output <- tibble::tibble(
    Analysis = 1:length(N), 
    N = N,
    rd = d * (p_c - p_e),
    rd0 = if(theta0 > 0){d * (p_c0 - p_e0)}else{0},
    theta = rd / sqrt(sigma2_H1),
    theta0 = rd0 / sqrt(sigma2_H0),
    info = info1, 
    info0 = info0
  )
  return(output)
}
