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

#'
#' @importFrom mvtnorm GenzBretz
#' @section Specification:
#' \if{latex}{
#'  \itemize{
#'    \item Define the analysis time from input fh_test.
#'    \item Compute arm0 and arm1 using \code{gs_create_arm()}.
#'    \item Set a unique test.
#'    \item Compute the information fraction using \code{gs_info_combo()}.
#'    \item Compute the correlation between tests.
#'    \item Compute the correlation between analysis.
#'    \item Compute the overall correlation.
#'    \item Extract the sample size from  info.
#'    \item Compute information restricted to actual analysis.
#'    \item Compute the effect size.
#'    \item Return a list of info_all = info, info = info_fh, theta = theta_fh, corr = corr_fh.
#'   }
#' }
#' \if{html}{The contents of this section are shown in PDF user manual only.}
#'
gs_utility_combo <- function(enrollRates,
                             failRates,
                             fh_test,
                             ratio = 1,
                             algorithm = GenzBretz(maxpts= 1e5, abseps= 1e-5),
                             ...){

  # Define analysis time
  analysisTimes <- sort(unique(fh_test$analysisTimes))

  # Define Arm
  gs_arm <- gs_create_arm(enrollRates, failRates,
                          ratio = ratio,                   # Randomization ratio
                          total_time = max(analysisTimes)) # Total study duration

  arm0 <- gs_arm[["arm0"]]
  arm1 <- gs_arm[["arm1"]]

  # Unique test
  u_fh_test <- unique(fh_test[, c("test","rho", "gamma", "tau")] )

  # Information Fraction
  info <- gs_info_combo(enrollRates, failRates, ratio,
                        analysisTimes = analysisTimes,
                        rho = u_fh_test$rho,
                        gamma = u_fh_test$gamma)

  # Correlation between test
  corr_test <- with(u_fh_test,
                    lapply(analysisTimes, function(tmax){
                      cov2cor(gs_sigma2_combo(arm0, arm1, tmax = tmax,
                                              rho = rho, gamma = gamma, tau = tau))
                    })
  )

  # Correlation between analysis
  info_split <- split(info, info$test)
  corr_time <- lapply(info_split, function(x){
    corr <- with(x, outer(sqrt(info), sqrt(info), function(x,y) pmin(x,y) / pmax(x,y)))
    rownames(corr) <- analysisTimes
    colnames(corr) <- analysisTimes
    corr
  })

  # Overall Correlation
  corr_combo <- diag(1, nrow = nrow(info))
  for(i in 1:nrow(info)){
    for(j in 1:nrow(info)){
      t1 <- as.numeric(info$Analysis[i])
      t2 <- as.numeric(info$Analysis[j])
      if(t1 <= t2){
        test1 <- as.numeric(info$test[i])
        test2 <- as.numeric(info$test[j])
        corr_combo[i,j] <- corr_test[[t1]][test1,test2] * corr_time[[test2]][t1, t2]
        corr_combo[j,i] <- corr_combo[i,j]
      }
    }
  }

  # Sample size
  n <- max(info$N)

  # Restricted to actual analysis
  info_fh <- merge(info, fh_test, all = TRUE)
  corr_fh <- corr_combo[! is.na(info_fh$gamma), ! is.na(info_fh$gamma)]
  info_fh <- subset(info_fh, ! is.na(gamma))

  # Effect size
  theta_fh <- (- info_fh$delta) / sqrt(info_fh$sigma2)

  list(info_all = info, info = info_fh, theta = theta_fh, corr = corr_fh)

}
