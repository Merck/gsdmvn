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

#' Group sequential design power under risk difference
#'
#' @param p_c rate at the control group
#' @param p_e rate at the experimental group 
#' @param n sample size 
#' @param theta0 the standardized treatment effect under H0  
#' @param delta0 treatment effect under super-superiority designs, the default is 0
#' @param ratio experimental:control randomization ratio
#' @param upper function to compute upper bound
#' @param upar parameter to pass to upper
#' @param lower function to compare lower bound
#' @param lpar parameter to pass to lower
#' @param info_scale 
#' @param binding indicator of whether futility bound is binding; default of FALSE is recommended
#' @param test_upper indicator of which analyses should include an upper (efficacy) bound;
#' single value of TRUE (default)  indicates all analyses; otherwise,
#' a logical vector of the same length as \code{info} should indicate which analyses will have an efficacy bound
#' @param test_lower indicator of which analyses should include a lower bound;
#' single value of TRUE (default) indicates all analyses;
#' single value FALSE indicated no lower bound; otherwise,
#' a logical vector of the same length as \code{info} should indicate which analyses will have a lower bound
#' @param r Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param tol Tolerance parameter for boundary convergence (on Z-scale)
#'
#' @return
#' @export
#'
#' @examples
#' gs_power_rd()
gs_power_rd <- function(
  p_c = .15,
  p_e = .13,
  N = c(50, 80, 100),
  theta0 = 0,
  delta0 = 0, 
  ratio = 1,
  upper = gs_b,
  upar = list(par = gsDesign(k = length(N), test.type = 1, sfu = sfLDOF, sfupar = NULL)$upper$bound),
  lower = gs_b,
  lpar = list(par = c(qnorm(.1), rep(-Inf, length(N) - 1))),
  info_scale = c(0, 1, 2),
  binding = FALSE,
  test_upper = TRUE,
  test_lower = TRUE,
  r = 18,
  tol = 1e-6
  ){
  
  # get the number of analysis
  K <- length(N)
  # get the info_scale
  info_scale <- if(methods::missingArg(info_scale)){2}else{match.arg(as.character(info_scale), choices = 0:2)}
  
  # ---------------------------------------- #
  #    calculate the asymptotic variance     #
  #       and statistical information        #
  # ---------------------------------------- #
  x <- gs_info_rd(
    p_c = p_c,
    p_e = p_e,
    N = N,
    theta0 = theta0,
    delta0 = delta0,
    ratio = ratio)
  
  # ---------------------------------------- #
  #  given the above statistical information #
  #         calculate the power              #
  # ---------------------------------------- #
  y_H1 <- gs_power_npe(
    theta = x$theta, 
    info = x$info, 
    info_scale = info_scale,
    binding = binding,
    upper = upper, 
    upar = upar,
    test_upper = test_upper,
    lower = lower, 
    lpar = lpar,
    test_lower = test_lower,
    r = r, 
    tol = tol) 
  
  y_H0 <- gs_power_npe(
    theta = x$theta0, 
    info = x$info0, 
    info_scale = info_scale,
    binding = binding,
    upper = upper, 
    upar = upar,
    test_upper = test_upper,
    lower = lower, 
    lpar = lpar,
    test_lower = test_lower,
    r = r, 
    tol = tol)
  
  # ---------------------------------------- #
  #         organize the outputs             #
  # ---------------------------------------- #
  # summarize the bounds
  suppressMessages(
    bounds <- y_H0 %>% 
      mutate(`~Risk difference at bound` = exp(-Z / sqrt(info)),  `Nominal p` = pnorm(-Z)) %>% 
      dplyr::rename(Probability0 = Probability) %>% 
      left_join(y_H1 %>% select(Analysis, Bound, Probability)) %>% 
      select(Analysis, Bound, Probability, Probability0, Z, `~Risk difference at bound`, `Nominal p`)
  )
  # summarize the analysis
  suppressMessages(
    analysis <- x %>% 
      select(Analysis, N, rd, rd0, theta, theta0) %>% 
      left_join(y_H1 %>% select(Analysis, info, IF) %>% unique()) %>%
      left_join(y_H0 %>% select(Analysis, info, IF) %>% dplyr::rename(info0 = info, IF0 = IF) %>% unique()) %>%
      select(Analysis, N, rd, rd0, theta, theta0, info, info0, IF, IF0)
  )
  
  output <- list(
    bounds = bounds,
    analysis = analysis)
  
  class(output) <- c("rd", "gs_design", class(output))
  
  return(output)
}