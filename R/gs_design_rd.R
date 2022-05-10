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
  N = tibble::tibble(Stratum = "All",  
                     N = c(20, 30, 50),
                     Analysis = 1:3),
  rd0 = 0, 
  alpha = 0.025,                   # One-sided Type I error
  beta = 0.1,                      # NULL if enrollment is not adapted
  ratio = 1,
  weight = c("un-stratified", "ss", "invar"),
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
  K <- length(N$N)
  info_scale <- if(methods::missingArg(info_scale)){2}else{match.arg(as.character(info_scale), choices = 0:2)}
  weight <- if(methods::missingArg(weight)){"un-stratified"}else{match.arg(weight)}
  
  # --------------------------------------------- #
  #     calculate the sample size                 #
  #          under fixed design                   #
  # --------------------------------------------- #
  N_fixed_des <- gsDesign::nBinomial(p1 = p_c$Rate, p2 = p_e$Rate, alpha = alpha, beta = beta, ratio = ratio, delta0 = rd0)
  
  # ---------------------------------------- #
  #    calculate the asymptotic variance     #
  #       and statistical information        #
  # ---------------------------------------- #
  x <- gs_info_rd(
    p_c = p_c,
    p_e = p_e,
    N = tibble::tibble(Stratum = "All",  
                       N = N_fixed_des *  N$N[1:length(N$N)] / N$N[length(N$N)],
                       Analysis = 1:K),
    rd0 = rd0,
    ratio = ratio,
    weight = weight)
  
  # --------------------------------------------- #
  #     get statistical information               #
  # --------------------------------------------- #
  suppressMessages(
    allout <- gs_design_npe(
      theta = x$theta, 
      info = x$info, info0 = x$info0, info_scale = info_scale,
      alpha = alpha, beta = beta, binding = binding,
      upper = upper, upar = upar, test_upper = test_upper,
      lower = lower, lpar = lpar, test_lower = test_lower,
      r = r, tol = tol) %>%
      # add `~HR at bound`, `HR generic` and `Nominal p`
      mutate("~Risk difference at bound" = exp(-Z / sqrt(info)), 
             "Nominal p" = pnorm(-Z),
             IF0 = if(sum(!is.na(info0)) == 0){NA}else{info0 / max(info0)}) %>% 
      # Add `Time`, `Events`, `AHR`, `N` from gs_info_ahr call above
      full_join(x %>% select(-c(info, info0, theta)), by = "Analysis") %>%
      # select variables to be output
      select(c(Analysis, Bound,  N, rd, rd0, Z, Probability, Probability0, theta, theta0, info, info0, IF, IF0, `~Risk difference at bound`, `Nominal p`)) %>% 
      # arrange the output table
      arrange(desc(Bound), Analysis) 
  )
  allout$N <- allout$N * allout$info[K] / x$info[K]
  
  # --------------------------------------------- #
  #     get bounds to output                      #
  # --------------------------------------------- #
  bounds <- allout %>%  
    select(Analysis, Bound, Probability, Probability0, Z, `~Risk difference at bound`, `Nominal p`)
  # --------------------------------------------- #
  #     get analysis summary to output            #
  # --------------------------------------------- #
  analysis <- allout %>% 
    select(Analysis, N, rd, rd0, theta, theta0, info, info0, IF, IF0) 
  
  # --------------------------------------------- #
  #     return the output                         #
  # --------------------------------------------- #
  output <- list(
    bounds = bounds,
    analysis = analysis)
  
  class(output) <- c("rd", "gs_design", class(output))
  
  return(output)
}
