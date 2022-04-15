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
#' 
#' @param p_c 
#' @param p_e 
#' @param n sample size 
#' @param delta 
#' @param delta0 
#' @param ratio Experimental:Control randomization ratio
#'
#' @export
#' @example gs_info_rd()
#' 
gs_info_rd <- function(
  p_c = .15,
  p_e = .1,
  n = c(350, 700, 1400),
  delta = NULL,
  delta0 = 0, 
  ratio = 1  
){
  
  n_e <- n / (1 + ratio)
  n_c <- n * ratio / (1 + ratio)
  n <- n_c + n_e
  p_pool <- p_c * n_c / n + p_e * n_e /n
  p_c_star <- p_c + delta * ratio / (1 + ratio)
  p_e_star <- p_c - delta * 1 / (1 + ratio)
  
  if (is.null(delta)) delta <- p_e - p_c
  theta <- rep(p_e - p_c, length(n))
  
  sigma2 <- p_c * (1 - p_c) / n_c  + p_e * (1 - p_e) / n_e
  sigma2_H0 <- p_pool * (1 - p_pool) * (1 / n_c + 1 / n_e)
  sigma2_H1 <- p_c_star * (1 - p_c_star) / n_c  + p_e_star * (1 - p_e_star) / n_e
  info <- 1 / sigma2 
  info0 <- 1 / sigma2_H0
  info1 <- 1 / sigma2_H1
  
  output <- tibble::tibble(
    Analysis = 1:length(n), n = n,
    theta = theta, theta1 = rep(delta, length(n)),
    info = info, info0 = info0, info1 = info1)
  return(output)
}
