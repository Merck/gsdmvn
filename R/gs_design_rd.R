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
#' @param rd0 treatment effect under super-superiority designs, the default is 0
#' @param info_scale 
#' @param weight the weighting scheme for stratified population
#' @param tol Tolerance parameter for boundary convergence (on Z-scale)
#'
#' @return a \code{tibble} with columns Analysis, Bound, Z, Probability, theta, Time, AHR, Events
#' @details Need to be added
#' @export
#'
#' @examples
#' gs_design_rd()
gs_design_rd <- function(
  p_c = tibble::tibble(Stratum = "All", 
                       Rate = .2),
  p_e = tibble::tibble(Stratum = "All",
                       Rate = .15),
  k = 3,
  IF = 1:k/k,
  rd0 = 0, 
  alpha = .025,                  
  beta = .1,                    
  ratio = 1,
  weight = c("un-stratified", "ss", "invar"),
  upper = gs_b,
  lower = gs_b,
  upar = list(par = gsDesign(k = k, test.type = 1, sfu = sfLDOF, sfupar = NULL)$upper$bound),
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
  weight <- if(methods::missingArg(weight)){"un-stratified"}else{match.arg(weight)}
  
  # --------------------------------------------- #
  #     calculate the sample size                 #
  #          under fixed design                   #
  # --------------------------------------------- #
  x_fix <- gs_info_rd(
    p_c = p_c, 
    p_e = p_e,
    N = tibble::tibble(Stratum = p_c$Stratum, N = 1, Analysis = 1),
    rd0 = rd0,
    ratio = ratio,
    weight = weight) 
  
  y_fix <- gs_design_npe(
    theta = x_fix$rd, 
    info = x_fix$info, 
    info0 = x_fix$info0, 
    info_scale = info_scale,
    alpha = alpha, beta = beta, binding = binding,
    upper = upper, upar = upar, test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol)
  
  # --------------------------------------------- #
  #     calculate the sample size                 #
  #     under group sequential design             #
  # --------------------------------------------- #
  x_gs <- gs_info_rd(
    p_c = p_c, 
    p_e = p_e,
    N = tibble::tibble(Stratum = p_c$Stratum, N = 1:k/k, Analysis = 1:k),
    rd0 = rd0,
    ratio = ratio,
    weight = weight)

  y_gs <- gs_design_npe(
    theta = x_gs$rd, 
    info = x_gs$info, 
    info0 = x_gs$info0, 
    info_scale = info_scale,
    alpha = alpha, beta = beta, binding = binding,
    upper = upper, upar = upar, test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol)
  
  # --------------------------------------------- #
  #     get statistical information               #
  # --------------------------------------------- #
  allout <-  y_gs%>%
    mutate(rd = x_fix$rd,
           rd0 = rd0,
           "~Risk difference at bound" = Z / sqrt(info) / theta * (rd -rd0)  + rd0, 
           "Nominal p" = pnorm(-Z),
           IF0 = if(sum(!is.na(info0)) == 0){NA}else{info0 / max(info0)},
           N = y_gs$info[k] / x_fix$info[1]  * IF) %>% 
    select(c(Analysis, Bound,  N, rd, rd0, Z, Probability, Probability0, info, info0, IF, IF0, `~Risk difference at bound`, `Nominal p`)) %>% 
    arrange(desc(Bound), Analysis) 
  
  # --------------------------------------------- #
  #     get bounds to output                      #
  # --------------------------------------------- #
  bounds <- allout %>%  
    select(Analysis, Bound, Probability, Probability0, Z, `~Risk difference at bound`, `Nominal p`)
  # --------------------------------------------- #
  #     get analysis summary to output            #
  # --------------------------------------------- #
  analysis <- allout %>% 
    select(Analysis, N, rd, rd0, info, info0, IF, IF0) 
  
  # --------------------------------------------- #
  #     return the output                         #
  # --------------------------------------------- #
  output <- list(
    bounds = bounds,
    analysis = analysis)
  
  class(output) <- c("rd", "gs_design", class(output))
  
  return(output)
}
