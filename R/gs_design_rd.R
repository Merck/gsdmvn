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
#' @importFrom gsDesign gsDesign sfLDOF
#' @importFrom stats qnorm
#' @importFrom gsDesign2 AHR
#' @importFrom dplyr mutate full_join select arrange desc
NULL
#' Group sequential design using average hazard ratio under non-proportional hazards
#'
#' @param ratio Experimental:Control randomization ratio (not yet implemented)
#' @param alpha One-sided Type I error
#' @param beta Type II error
#' @param binding indicator of whether futility bound is binding; default of FALSE is recommended
#' @param upper Function to compute upper bound
#' @param upar Parameter passed to \code{upper()}
#' @param lower Function to compute lower bound
#' @param lpar Parameter passed to \code{lower()}
#' @param test_upper indicator of which analyses should include an upper (efficacy) bound; single value of TRUE (default) indicates all analyses;
#' otherwise, a logical vector of the same length as \code{info} should indicate which analyses will have an efficacy bound
#' @param test_lower indicator of which analyses should include an lower bound; single value of TRUE (default) indicates all analyses;
#' single value FALSE indicated no lower bound; otherwise, a logical vector of the same length as \code{info} should indicate which analyses will have a
#' lower bound
#' @param r  Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param p_c rate at the control group
#' @param p_e rate at the experimental group 
#' @param n sample size 
#' @param theta0 the standardized treatment effect under H0 
#' @param delta0 treatment effect under super-superiority designs, the default is 0
#' @param info_scale 
#' @param tol Tolerance parameter for boundary convergence (on Z-scale)
#'
#' @return a \code{tibble} with columns Analysis, Bound, Z, Probability, theta, Time, AHR, Events
#' @details Need to be added
#' @export
#'
#' @examples
#' gs_design_rd()
gs_design_rd <- function(
  p_c = .15,
  p_e = .13,
  N = c(2, 4, 6),
  theta0 = 0,
  delta0 = 0, 
  alpha = 0.025,                   # One-sided Type I error
  beta = 0.1,                      # NULL if enrollment is not adapted
  ratio = 1,
  upper = gs_b,
  lower = gs_b,
  upar = list(par = gsDesign(k = length(N), test.type = 1, sfu = sfLDOF, sfupar = NULL)$upper$bound),
  lpar = list(par = c(qnorm(.1), rep(-Inf, length(N) - 1))),
  test_upper = TRUE,
  test_lower = TRUE,
  info_scale = c(0, 1, 2),
  binding = FALSE,
  r = 18,
  tol = 1e-6
){
  # --------------------------------------------- #
  #     check input values                        #
  # --------------------------------------------- #
  info_scale <- if(methods::missingArg(info_scale)){2}else{match.arg(as.character(info_scale), choices = 0:2)}
  K <- length(N)
  # --------------------------------------------- #
  #     get statistical information               #
  # --------------------------------------------- #
  x <- gs_info_rd(
    p_c = p_c,
    p_e = p_e,
    N = N,
    theta0 = theta0,
    delta0 = delta0,
    ratio = ratio)
  
  find_xi <- function(xi){
    x_temp <- gs_power_rd(
      p_c = p_c,
      p_e = p_e,
      N = N * xi,
      theta0 = theta0,
      delta0 = delta0, 
      ratio = ratio,
      upper = upper,
      upar = upar,
      lower = lower,
      lpar = lpar,
      info_scale = info_scale,
      binding = binding,
      test_upper = test_upper,
      test_lower = test_lower,
      r = r,
      tol = tol)
    
    x_temp_power <- x_temp$bounds %>% 
      filter(Bound == "Upper", Analysis == K) %>% 
      select(Probability) %>% 
      unlist() %>% 
      as.numeric()
    
    return(1 - beta - x_temp_power)
  }
  
  res <- try(uniroot(find_xi, lower = 1, upper = 50))
  if(inherits(res, "try-error")){stop("gs_design_rd: Sample size solution not found")}
  if(abs(res$f.root) >  1e-5){stop("gs_design_rd: Sample size solution not found")}
  # --------------------------------------------- #
  #     combine all the calculations              #
  # --------------------------------------------- #
  x_final <- gs_power_rd(
    p_c = p_c, p_e = p_e,
    N = N * res$root,
    theta0 = theta0,
    delta0 = delta0, 
    ratio = ratio,
    upper = upper,
    upar = upar,
    lower = lower,
    lpar = lpar,
    info_scale = info_scale,
    binding = binding,
    test_upper = test_upper,
    test_lower = test_lower,
    r = r,
    tol = tol
  )
  
  # --------------------------------------------- #
  #     get bounds to output                      #
  # --------------------------------------------- #
  bounds <- x_final$bounds %>% arrange(desc(Bound), Analysis)
  # --------------------------------------------- #
  #     get analysis summary to output            #
  # --------------------------------------------- #
  analysis <- x_final$analysis %>% arrange(Analysis)
  
  # --------------------------------------------- #
  #     return the output                         #
  # --------------------------------------------- #
  output <- list(
    bounds = bounds,
    analysis = analysis)
  class(output) <- c("rd", "gs_design", class(output))
  return(output)
}
