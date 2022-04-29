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
#' @param delta0 the treatment effect under H0 
#' @param ratio Experimental:Control randomization ratio
#'
#' @export
#' @examples
#' # --------------------- #
#' #      example 1        #
#' # --------------------- #
#' # un-stratified case with H0: theta = 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = "All",
#'                        Rate = c(.15, .18, .2),
#'                        Analysis = 1:3),
#'   p_e = tibble::tibble(Stratum = "All",
#'                        Rate = c(.12, .13, .15),
#'                        Analysis = 1:3),
#'   N = tibble::tibble(Stratum = "All",
#'                      N = c(100, 200, 300),
#'                      Analysis = 1:3),
#'   delta0 = 0,
#'   ratio = 1
#' )
#' 
#' # --------------------- #
#' #      example 2        #
#' # --------------------- #
#' # un-stratified case with H0: theta != 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = "All",
#'                        Rate = c(.15, .18, .2),
#'                        Analysis = 1:3),
#'   p_e = tibble::tibble(Stratum = "All",
#'                        Rate = c(.12, .13, .15),
#'                        Analysis = 1:3),
#'   N = tibble::tibble(Stratum = "All",
#'                      N = c(100, 200, 300),
#'                      Analysis = 1:3),
#'   delta0 = 0.005,
#'   ratio = 1
#' )
#' 
#' # --------------------- #
#' #      example 3        #
#' # --------------------- #
#' # stratified case under sample size weighting and H0: theta = 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.15, .2, .25, .13, .15, .18, .17, .22, .28)),
#'   p_e = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.12, .18, .21, .1, .11, .12, .11, .19, .26)),
#'   N = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                      Analysis = rep(1:3, 3),
#'                      N = c(50, 100, 200, 40, 80, 160, 60, 120, 240)),
#'   delta0 = 0,
#'   ratio = 1,
#'   weight = "ss")
#' 
#' # --------------------- #
#' #      example 4        #
#' # --------------------- #
#' # stratified case under inverse variance weighting and H0: theta = 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.15, .2, .25, .13, .15, .18, .17, .22, .28)),
#'   p_e = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.12, .18, .21, .1, .11, .12, .11, .19, .26)),
#'   N = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                      Analysis = rep(1:3, 3),
#'                      N = c(50, 100, 200, 40, 80, 160, 60, 120, 240)),
#'   delta0 = 0,
#'   ratio = 1,
#'   weight = "invar")
#' 
#' # --------------------- #
#' #      example 5        #
#' # --------------------- #
#' # stratified case under sample size weighting and H0: theta != 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.15, .2, .25, .13, .15, .18, .17, .22, .28)),
#'   p_e = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.12, .18, .21, .1, .11, .12, .11, .19, .26)),
#'   N = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                      Analysis = rep(1:3, 3),
#'                      N = c(50, 100, 200, 40, 80, 160, 60, 120, 240)),
#'   delta0 = 0,
#'   ratio = 1,
#'   weight = "ss")
#' 
#' # --------------------- #
#' #      example 6        #
#' # --------------------- #
#' # stratified case under inverse variance weighting and H0: theta != 0
#' gs_info_rd(
#'   p_c = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.15, .2, .25, .13, .15, .18, .17, .22, .28)),
#'   p_e = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                        Analysis = rep(1:3, 3),
#'                        Rate = c(.12, .18, .21, .1, .11, .12, .11, .19, .26)),
#'   N = tibble::tibble(Stratum = rep(c("S1", "S2", "S3"), each = 3),
#'                      Analysis = rep(1:3, 3),
#'                      N = c(50, 100, 200, 40, 80, 160, 60, 120, 240)),
#'   delta0 = 0,
#'   ratio = 1,
#'   weight = "invar")
#'   
gs_info_rd <- function(
  p_c = tibble::tibble(Stratum = rep(c("Strata 1", "Strata 2", "Strata 3"), each = 3), 
                       Analysis = rep(1:3, 3), 
                       Rate = c(.15, .2, .25, .13, .15, .18, .17, .22, .28)),
  p_e = tibble::tibble(Stratum = rep(c("Strata 1", "Strata 2", "Strata 3"), each = 3), 
                       Analysis = rep(1:3, 3), 
                       Rate = c(.12, .18, .21, .1, .11, .12, .11, .19, .26)),
  N = tibble::tibble(Stratum = rep(c("Strata 1", "Strata 2", "Strata 3"), each = 3),  
                     Analysis = rep(1:3, 3), 
                     N = c(50, 100, 200, 40, 80, 160, 60, 120, 240)),
  delta0 = 0, 
  ratio = 1,
  weight = c("ss", "invar")
){
  
  K <- max(N$Analysis)
  weight <- if(methods::missingArg(weight)){"un-stratified"}else{match.arg(weight)}
  
  # pool the input arguments together
  suppressMessages(
    tbl <- N %>% 
      left_join(p_c) %>%
      dplyr::rename(p_c = Rate) %>% 
      left_join(p_e) %>% 
      dplyr::rename(p_e = Rate) %>% 
      mutate(
        N_e = N / (1 + ratio), N_c = N * ratio / (1 + ratio),
        n_c = N_c * p_c, n_e = N_e * p_e, n = n_c + n_e,
        d = ifelse(p_c > p_e, 1, -1),
        p_pool_per_stratum = (n_c + n_e) / N,
        p_e0 = (p_c + ratio * p_e - d * delta0) / (ratio + 1),
        p_c0 = p_e0 + d * delta0)
  )
  
  # calculate the variance of the risk difference
  if(delta0 == 0){
    tbl <- tbl %>% mutate(sigma2_H0 = p_pool_per_stratum * (1 - p_pool_per_stratum) * (1 / N_c + 1 / N_e),
                          sigma2_H1 = p_c * (1 - p_c) / N_c  + p_e * (1 - p_e) / N_e)  
  }
  if(delta0 > 0){
    tbl <- tbl %>% mutate(sigma2_H0 = p_c0 * (1 - p_c0) / N_c + p_e0 * (1 - p_e0) / N_e,
                          sigma2_H1 = p_c * (1 - p_c) / N_c  + p_e * (1 - p_e) / N_e) 
  }

  # assign weights
  if(weight == "un-stratified"){
    tbl <- tbl %>% mutate(weight = 1) 
  }else if(weight == "ss"){
    suppressMessages(
      tbl <- tbl %>% 
        left_join(tbl %>% dplyr::group_by(Analysis) %>% summarize(sum_ss = sum(N_c * N_e / (N_c + N_e)))) %>%
        mutate(weight = N_c * N_e / (N_c + N_e) / sum_ss ) %>% 
        select(-sum_ss)
    )
  }else if(weight == "invar"){
    suppressMessages(
      tbl <- tbl %>% 
        left_join(tbl %>% dplyr::group_by(Analysis) %>% summarize(sum_inv_var = sum(1/sigma2_H0))) %>% 
        mutate(weight = 1/sigma2_H0 / sum_inv_var) %>% 
        select(-sum_inv_var)
    )
  }
  
  # pool risk difference, theta, sigma together from different stratum together
  if(delta0 == 0){
    tbl <- tbl %>% mutate(sigma2_H0_per_strata = weight^2 * p_pool_per_stratum * (1 - p_pool_per_stratum) / (N_c + N_e),
                          sigma2_H1_per_strata = weight^2 * p_c * (1 - p_c) / N_c + weight^2 * p_e * (1 - p_e) / N_e) 
  }
  
  if(delta0 > 0){
    tbl <- tbl %>% mutate(sigma2_H0_per_strata = p_c0 * (1 - p_c0) / N_c +  p_e0 * (1 - p_e0) / N_e,
                          sigma2_H1_per_strata = p_c * (1 - p_c) / N_c  + p_e * (1 - p_e) / N_e)
  }

  output <- tbl %>% 
    group_by(Analysis) %>%
    summarize(rd = sum((p_c - p_e) * d * weight),
              rd0 = sum(delta0 * weight),
              N = sum(N), N_c = sum(N_c), N_e = sum(N_e),
              n = sum(n), n_c = sum(n_c), n_e = sum(n_e),
              sigma2_H0 = sum(weight * sigma2_H0_per_strata),
              sigma2_H1 = sum(weight * sigma2_H1_per_strata),
              theta = rd / sqrt(sigma2_H1),
              theta0 = rd0 / sqrt(sigma2_H0),
              info = 1 / sigma2_H1,
              info0 = 1 / sigma2_H0) %>%
    ungroup() %>% 
    select(Analysis, N, rd, rd0, theta, theta0, info, info0)

  return(output)
}
