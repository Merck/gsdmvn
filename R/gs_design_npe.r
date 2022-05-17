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
#' @importFrom stats qnorm uniroot
NULL
#' Group sequential design computation with non-constant effect and information
#'
#' \code{gs_design_npe()} derives group sequential design size, bounds and boundary crossing probabilities based on proportionate
#' information and effect size at analyses.
#' It allows a non-constant treatment effect over time, but also can be applied for the usual homogeneous effect size designs.
#' It requires treatment effect and proportionate statistical information at each analysis as well as a method of deriving bounds, such as spending.
#' The routine enables two things not available in the gsDesign package: 1) non-constant effect, 2) more flexibility in boundary selection.
#' For many applications, the non-proportional-hazards design function \code{gs_design_nph()} will be used; it calls this function.
#' Initial bound types supported are 1) spending bounds, 2) fixed bounds, and 3) Haybittle-Peto-like bounds.
#' The requirement is to have a boundary update method that can each bound without knowledge of future bounds.
#' As an example, bounds based on conditional power that require knowledge of all future bounds are not supported by this routine;
#' a more limited conditional power method will be demonstrated.
#' Boundary family designs Wang-Tsiatis designs including the original (non-spending-function-based) O'Brien-Fleming and Pocock designs
#' are not supported by \code{gs_power_npe()}.
#' @param theta natural parameter for group sequential design representing expected incremental drift at all analyses;
#' used for power calculation
#' @param theta1 natural parameter used for lower bound spending; if \code{NULL}, this will be set to \code{theta}
#' which yields the usual beta-spending. If set to 0, spending is 2-sided under null hypothesis.
#' @param info proportionate statistical information at all analyses for input \code{theta}
#' @param info0 proportionate statistical information under null hypothesis, if different than alternative;
#' impacts null hypothesis bound calculation
#' @param alpha One-sided Type I error
#' @param beta Type II error
#' @param binding indicator of whether futility bound is binding; default of FALSE is recommended
#' @param upper function to compute upper bound
#' @param lower function to compare lower bound
#' @param upar parameter to pass to function provided in \code{upper}
#' @param lpar Parameter passed to function provided in \code{lower}
#' @param test_upper indicator of which analyses should include an upper (efficacy) bound; single value of TRUE (default) indicates all analyses;
#' otherwise, a logical vector of the same length as \code{info} should indicate which analyses will have an efficacy bound
#' @param test_lower indicator of which analyses should include an lower bound; single value of TRUE (default) indicates all analyses;
#' single value FALSE indicated no lower bound; otherwise, a logical vector of the same length as \code{info} should indicate which analyses will have a
#' lower bound
#' @param r  Integer, at least 2; default of 18 recommended by Jennison and Turnbull
#' @param tol Tolerance parameter for boundary convergence (on Z-scale)
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Validate if input info is a numeric vector  or NULL, if non-NULL validate if it
#'    is strictly increasing and positive.
#'    \item Validate if input info0 is a numeric vector or NULL, if non-NULL validate if it
#'     is strictly increasing and positive.
#'    \item Validate if input info1 is a numeric vector or NULL, if non-NULL validate if it
#'    is strictly increasing and positive.
#'    \item Validate if input theta is a real vector and has the same length as info.
#'    \item Validate if input theta1 is a real vector and has the same length as info.
#'    \item Validate if input test_upper and test_lower are logical and have the same length as info.
#'    \item Validate if input test_upper value is TRUE.
#'    \item Validate if input alpha and beta are positive and of length one.
#'    \item Validate if input alpha and beta are from the unit interval and alpha is smaller than beta.
#'    \item Initialize bounds, numerical integration grids, boundary crossing probabilities.
#'    \item Compute fixed sample size for desired power and Type I error.
#'    \item Find an interval for information inflation to give correct power using \code{gs_power_npe()}.

#'    \item
#'    \item If there is no interim analysis, return a tibble including Analysis time, upper bound, Z-value,
#'    Probability of crossing bound, theta, info0 and info1.
#'    \item If the design is a group sequential design, return a tibble of Analysis,
#'     Bound, Z, Probability,  theta, info, info0.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
#' @return a \code{tibble} with columns Analysis, Bound, Z, Probability,  theta, info, info0
#' @details The inputs \code{info} and \code{info0} should be vectors of the same length with increasing positive numbers.
#' The design returned will change these by some constant scale factor to ensure the design has power \code{1 - beta}.
#' The bound specifications in \code{upper, lower, upar, lpar} will be used to ensure Type I error and other boundary properties are as specified.
#' @author Keaven Anderson \email{keaven_anderson@@merck.com}
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(gsDesign)
#' 
#' # ---------------------------------# 
#' #         example 1                #
#' # ---------------------------------# 
#' # Single analysis
#' # Lachin book p 71 difference of proportions example
#' pc <- .28            # Control response rate
#' pe <- .40            # Experimental response rate
#' p0 <- (pc + pe) / 2  # Ave response rate under H0
#' 
#' # Information per increment of 1 in sample size
#' info0 <- 1 / (p0 * (1 - p0) * 4)
#' info <- 1 / (pc * (1 - pc) * 2 + pe * (1 - pe) * 2)
#' 
#' # Result should round up to next even number = 652
#' # Divide information needed under H1 by information per patient added
#' gs_design_npe(theta = pe - pc, info = info, info0 = info0)
#' 
#' # One can try `info_scale` argument. But it gives the same results as above. 
#' # This is because the above example use fixed design.
#' gs_design_npe(theta = pe - pc, info = info, info0 = info0, info_scale = 0)
#' gs_design_npe(theta = pe - pc, info = info, info0 = info0, info_scale = 1)
#' gs_design_npe(theta = pe - pc, info = info, info0 = info0, info_scale = 2) # default 
#' 
#' # ---------------------------------# 
#' #         example 2                #
#' # ---------------------------------# 
#' # Fixed bound
#' x <- gs_design_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 80,
#'   info0 = (1:3) * 80,
#'   upper = gs_b,
#'   upar = list(par = gsDesign::gsDesign(k = 3, sfu = gsDesign::sfLDOF)$upper$bound),
#'   lower = gs_b,
#'   lpar = list(par = c(-1, 0, 0)))
#' x
#' 
#' # Same upper bound; this represents non-binding Type I error and will total 0.025
#' gs_power_npe(
#'   theta = rep(0, 3),
#'   info = (x %>% filter(Bound == "Upper"))$info,
#'   upper = gs_b,
#'   upar = list(par = (x %>% filter(Bound == "Upper"))$Z),
#'   lower = gs_b,
#'   lpar = list(par = rep(-Inf, 3)))
#' 
#' # ---------------------------------# 
#' #         example 3                #
#' # ---------------------------------# 
#' # Spending bound examples
#' # Design with futility only at analysis 1; efficacy only at analyses 2, 3
#' # Spending bound for efficacy; fixed bound for futility
#' # NOTE: test_upper and test_lower DO NOT WORK with gs_b; must explicitly make bounds infinite
#' # test_upper and test_lower DO WORK with gs_spending_bound
#' gs_design_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 40,
#'   info0 = (1:3) * 40,
#'   upper = gs_spending_bound,
#'   upar = list(par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)),
#'   lower = gs_b,
#'   lpar = list(par = c(-1, -Inf, -Inf)),
#'   test_upper = c(FALSE, TRUE, TRUE))
#' # one can try `info_scale = 1` or `info_scale = 0` here
#' gs_design_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 40,
#'   info0 = (1:3) * 30,
#'   info_scale = 1,
#'   upper = gs_spending_bound,
#'   upar = list(par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)),
#'   lower = gs_b,
#'   lpar = list(par = c(-1, -Inf, -Inf)),
#'   test_upper = c(FALSE, TRUE, TRUE))
#' 
#' # ---------------------------------# 
#' #         example 4                #
#' # ---------------------------------# 
#' # Spending function bounds
#' # 2-sided asymmetric bounds
#' # Lower spending based on non-zero effect
#' gs_design_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 40,
#'   info0 = (1:3) * 30,
#'   upper = gs_spending_bound,
#'   upar = list(par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL),
#'               info = (1:3) * 30),
#'   lower = gs_spending_bound,
#'   lpar = list(par = list(sf = gsDesign::sfHSD, total_spend = 0.1, param = -1, timing = NULL)))
#' 
#' # ---------------------------------# 
#' #         example 5                #
#' # ---------------------------------# 
#' # Two-sided symmetric spend, O'Brien-Fleming spending
#' # Typically, 2-sided bounds are binding
#' xx <- gs_design_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 40,
#'   binding = TRUE,
#'   upper = gs_spending_bound,
#'   upar = list(par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)),
#'   lower = gs_spending_bound,
#'   lpar = list(par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL)))
#' xx
#' 
#' # Re-use these bounds under alternate hypothesis
#' # Always use binding = TRUE for power calculations
#' gs_power_npe(
#'   theta = c(.1, .2, .3),
#'   info = (1:3) * 40,
#'   binding = TRUE,
#'   upper = gs_b,
#'   lower = gs_b,
#'   upar = list(par = (xx %>% filter(Bound == "Upper"))$Z),
#'   lpar = list(par = -(xx %>% filter(Bound == "Upper"))$Z))
#'   
gs_design_npe <- function(
  theta = .1, 
  info = 1, 
  info0 = NULL, 
  info_scale = c(0, 1, 2),
  alpha = 0.025, beta = .1, 
  binding = FALSE,
  upper = gs_b, upar = list(par = qnorm(.975)), test_upper = TRUE,
  lower = gs_b, lpar = list(par = -Inf), test_lower = TRUE,
  r = 18, tol = 1e-6){
  
  # --------------------------------------------- #
  #     check info, info0                         #
  # --------------------------------------------- #
  if (!is.vector(info, mode = "numeric")) stop("gs_design_npe(): info must be specified numeric vector")
  K <- length(info)
  if (is.null(info0)) info0 <- info
  if (!is.vector(info0, mode = "numeric")) stop("gs_design_npe(): info0 must be specified numeric vector or NULL")
  if (length(info0) != length(info) ) stop("gs_design_npe(): length of info, info0 must be the same")
  if (min(info - lag(info, default = 0) <= 0)) stop("gs_design_npe(): info much be strictly increasing and positive")
  if (min(info0 - lag(info0, default = 0) <= 0)) stop("gs_design_npe(): info0 much be NULL or strictly increasing and positive")
  
  if(identical(upper, gs_spending_bound)){
    info_scale <- if(methods::missingArg(info_scale)){2}else{match.arg(as.character(info_scale), choices = 0:2)}
  }
  # --------------------------------------------- #
  #     check theta, theta0, theta1               #
  # --------------------------------------------- #
  if (!is.vector(theta, mode = "numeric")) stop("gs_design_npe(): theta must be a real vector")
  if (length(theta) == 1 && K > 1) theta <- rep(theta, K)
  if (length(theta) != K) stop("gs_design_npe(): if length(theta) > 1, must be same as info")
  if (theta[K] <= 0) stop("gs_design_npe(): final effect size must be > 0")
  
  # --------------------------------------------- #
  #     check test_upper & test_lower             #
  # --------------------------------------------- #
  # check the correct spec of test_upper & test_lower
  if (length(test_upper) == 1 && K > 1) test_upper <- rep(test_upper, K)
  if (length(test_lower) == 1 && K > 1) test_lower <- rep(test_lower, K)
  ## Check test_upper and test_lower are logical and correct length
  if (!is.vector(test_upper, mode = "logical") || !is.vector(test_lower, mode = "logical"))
    stop("gs_design_npe(): test_upper and test_lower must be logical")
  if (!(length(test_upper) == 1 || length(test_upper) == K))
    stop("gs_design_npe(): test_upper must be length 1 or same length as info")
  if (!(length(test_lower) == 1 || length(test_lower) == K))
    stop("gs_design_npe(): test_lower must be length 1 or same length as info")
  ## Check that final test_upper value is TRUE
  if (!dplyr::last(test_upper)) stop("gs_design_npe(): last value of test_upper must be TRUE")
  
  # --------------------------------------------- #
  #     check alpha & beta                        #
  # --------------------------------------------- #
  ## Check alpha and beta numeric, scalar, 0 < alpha < 1 - beta
  if (!is.numeric(alpha)) stop("gs_design_npe(): alpha must be numeric")
  if (!is.numeric(beta)) stop("gs_design_npe(): beta must be numeric")
  if (length(alpha) != 1 || length(beta) != 1) stop("gs_design_npe(): alpha and beta must be length 1")
  if (alpha <= 0 || 1 - beta <= alpha || beta <= 0) stop("gs_design_npe(): must have 0 < alpha < 1 - beta < 1")
  
  # --------------------------------------------- #
  #     initialization                            #
  # --------------------------------------------- #
  a <- rep(-Inf, K)          # bounds
  b <- rep(Inf, K)
  hgm1_0 <- NULL             # numerical integration grids
  hgm1_1 <- NULL
  upperProb <- rep(NA, K)    # boundary crossing probabilities
  lowerProb <- rep(NA, K)
  
  # --------------------------------------------- #
  #     fixed design                              #
  # --------------------------------------------- #
  # compute fixed sample size for desired power and Type I error.
  minx <- ((qnorm(alpha) / sqrt(info0[K]) + qnorm(beta) / sqrt(info[K])) / theta[K])^2
  # for a fixed design, this is all you need.
  if (K == 1){
    out <- tibble::tibble(Analysis = 1, Bound = "Upper", Z = qnorm(1 - alpha),
                          Probability = 1 - beta, Probability0 = alpha, theta = theta, 
                          info = info * minx, info0 = info0 * minx, IF = info / max(info))
    return(out)
  } 
  
  # find an interval for information inflation to give correct power
  minpwr <- gs_power_npe(
    theta = theta, 
    info = info * minx, 
    info_scale = info_scale,
    binding = binding,
    upper = upper, upar = c(upar, info = list(info0 * minx)), test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol
  ) %>% 
    filter(Bound == "Upper" & Analysis == K) %>% 
    select(Probability) %>% 
    unlist() %>% 
    as.numeric()  
  
  # --------------------------------------------- #
  #     FOLLOWING IS PAINFUL                      #
  #       AND SHOULD NEVER BE NEEDED              #
  #     BUT IF IT IS NEEDED,                      #
  #       IT TELLS YOU WHAT WENT WRONG!           #
  #     NEED TO BRACKET TARGETED POWER            #
  #       BEFORE ROOT FINDING                     #
  # --------------------------------------------- #
  ## Ensure minx gives power < 1 - beta and maxx gives power > 1 - beta
  if (minpwr < 1 - beta){
    ## Ensure maxx is sufficient information inflation to overpower
    maxx <- 1.05 * minx
    err <- 1
    for(i in 1:10){
      maxpwr <- gs_power_npe(
        theta = theta, 
        info = info * maxx, 
        info_scale = info_scale,
        binding = binding,
        upper = upper, upar = c(upar, info = list(info * maxx)), test_upper = test_upper, 
        lower = lower, lpar = lpar, test_lower = test_lower,
        r = r, tol = tol
      )%>% 
        filter(Bound == "Upper" & Analysis == K) %>% 
        select(Probability) %>% 
        unlist() %>% 
        as.numeric() 
      
      if (1  - beta > maxpwr){
        minx <- maxx
        maxx <- 1.05 * maxx
      }else{
        err <- 0
        break
      }
    }
    if (err) stop("gs_design_npe: could not inflate information to bracket power before root finding")
  }else{
    maxx <- minx
    minx <- maxx / 1.05
    err <- 1
    for(i in 1:10){
      if (1  - beta < gs_power_npe(
        theta = theta, 
        info = info * minx, 
        info_scale = info_scale,
        binding = binding,
        upper = upper, lower = lower, 
        upar = c(upar, info = list(info0 * minx)), lpar = lpar,
        test_upper = test_upper, test_lower = test_lower,
        r = r, tol = tol
      ) %>% 
      filter(Bound == "Upper" & Analysis == K) %>% 
      select(Probability) %>% 
      unlist() %>% 
      as.numeric() 
      ){
        maxx <- minx
        minx <- minx / 1.05
      }else{
        err <- 0
        break
      }
    }
    if (err) stop("gs_design_npe: could not deflate information to bracket targeted power before root finding")
  }
  
  # --------------------------------------------- #
  #     EITHER TARGETED POWER NOW BRACKETED       #
  #                   OR                          #
  #     ERROR MESSAGE HAS BEEN RETURNED           #
  #     AND WE CAN ACTUALLY GO ON TO FIND THE ROOT#
  # --------------------------------------------- #
  ## Use root finding with the above function to find needed sample size inflation
  # Now we can solve for the inflation factor for the enrollment rate to achieve the desired power
  res <- try(
    uniroot(errbeta, lower = minx, upper = maxx,
            theta = theta, 
            K = K, 
            beta = beta,
            info = info, info0 = info0, info_scale = info_scale,
            binding = binding,
            Zupper = upper, Zlower = lower, 
            upar = upar, lpar = lpar,
            test_upper = test_upper, test_lower = test_lower,
            r = r, tol = tol)
  )
  if(inherits(res, "try-error")){stop("gs_design_npe: Sample size solution not found")}
  
  
  # --------------------------------------------- #
  #     return the output                         #
  # --------------------------------------------- #
  # calculate the probability under H1
  out_H1 <- gs_power_npe(
    theta = theta, 
    info = info * res$root, 
    info_scale = info_scale,
    binding = binding,
    upper = upper, 
    lower = lower,
    upar = c(upar, info = list(info0 * res$root)), 
    lpar = lpar, 
    test_upper = test_upper, test_lower = test_lower,
    r = r, tol = tol)
  # get the bounds from out_H1
  suppressMessages(
    bound_H1 <- out_H1 %>% 
      select(Analysis, Bound, Z) %>%
      dplyr::rename(Z1 = Z) %>% 
      right_join(tibble::tibble(Analysis = rep(1:K, 2), Bound = rep(c("Upper", "Lower"), each = K), Z2 = rep(c(Inf, -Inf), each = K))) %>% 
      mutate(Z = dplyr::coalesce(Z1, Z2)) %>% 
      select(Analysis, Bound, Z) %>% 
      arrange(desc(Bound), Analysis)
  )
  # calculate the probability under H0
  out_H0 <- gs_power_npe(
    theta = 0, 
    info = info * res$root, 
    info_scale = info_scale,
    binding = binding,
    upper = gs_b, upar = list(par = (bound_H1 %>% filter(Bound == "Upper"))$Z), test_upper = test_upper,
    lower = gs_b, lpar = list(par = (bound_H1 %>% filter(Bound == "Lower"))$Z), test_lower = test_lower,
    r = r, tol = tol)
  # combine probability  under H0 and H1
  suppressMessages(
    out <- out_H1 %>% full_join(out_H0 %>% select(Analysis, Bound, Z, Probability) %>% dplyr::rename(Probability0 = Probability))
  )
  if ("info0" %in% colnames(out)){
    out <- out %>% select(Analysis, Bound, Z, Probability, Probability0, theta, IF, info, info0) 
  }else{
    out <- out %>% select(Analysis, Bound, Z, Probability, Probability0, theta, IF, info) 
  }
  out <- out %>% arrange(desc(Bound), Analysis)
  return(out)
  
  
}


## Create a function that uses gs_power_npe to compute difference from targeted power
## for a given sample size inflation factor
errbeta <- function(x = 1, K = 1, 
                    beta = .1, 
                    theta = .1, 
                    info = 1, 
                    info0 = 1,
                    info_scale = 2,
                    binding = FALSE,
                    Zupper = gs_b, upar = qnorm(.975), test_upper = TRUE,
                    Zlower = gs_b, lpar= -Inf, test_lower = TRUE,
                    r = 18, tol = 1e-6){
  out <- 1 -  
    beta -
    gs_power_npe(theta = theta, 
                 info = info * x, 
                 binding = binding,
                 info_scale = info_scale,
                 upper = Zupper, lower = Zlower, 
                 upar = c(upar, info = list(info0 * x)), lpar = lpar,
                 test_upper = test_upper, test_lower = test_lower,
                 r = r, tol = tol
    )%>% 
    filter(Bound == "Upper" & Analysis == K) %>% 
    select(Probability) %>% 
    unlist() %>% 
    as.numeric()
  
  return(out) 
}