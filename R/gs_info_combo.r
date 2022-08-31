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

#' @noRd
get_combo_weight <- function(rho, gamma, tau){

  stopifnot(length(rho) == length(gamma))
  stopifnot(length(rho) == length(tau))

  weight <- list()
  for(i in 1:length(rho)){

    if(tau[i] == -1){
      tmp_tau <- NULL
    }else{
      tmp_tau <- tau[i]
    }

    text <- paste0("weight <- function(x, arm0, arm1){
                            wlr_weight_fh(x, arm0, arm1
                                ,rho =", rho[i],
                   ", gamma =", gamma[i],
                   ", tau =", tmp_tau,  ")}")

    weight[[i]] <- text

  }

  weight
}

#' @noRd
gs_delta_combo <- function(arm0,
                           arm1,
                           tmax = NULL,
                           rho,
                           gamma,
                           tau = rep(-1, length(rho)),
                           approx="asymptotic",
                           normalization = FALSE) {

  stopifnot(length(tmax) == 1)

  weight <- get_combo_weight(rho, gamma, tau)
  delta <- sapply(weight, function(x){
    x <- eval(parse(text = x))
    gs_delta_wlr(arm0, arm1, tmax = tmax, weight = x,
                 approx = approx, normalization = normalization)
  })

  delta

}

#' @noRd
gs_sigma2_combo <- function(arm0,
                            arm1,
                            tmax = NULL,
                            rho,
                            gamma,
                            tau = rep(-1, length(rho)),
                            approx="asymptotic"){

  stopifnot(length(tmax) == 1)
  stopifnot(length(rho) == length(gamma))
  stopifnot(length(rho) == length(tau))

  rho1   <- outer(rho, rho, function(x,y) (x+y)/2 )
  gamma1 <- outer(gamma, gamma, function(x,y) (x+y)/2 )

  sigma2 <- rho1
  for(i in 1:length(rho)){

    weight <- get_combo_weight(rho1[i,], gamma1[i,],tau)

    sigma2[i,] <- sapply(weight, function(x){
                    x <- eval(parse(text = x))
                    gs_sigma2_wlr(arm0, arm1, tmax = tmax, weight = x,
                                  approx = approx)

                  })
  }

  sigma2

}

#' @noRd
gs_info_combo <- function(enrollRates=tibble::tibble(Stratum="All",
                                                     duration=c(2,2,10),
                                                     rate=c(3,6,9)),
                          failRates=tibble::tibble(Stratum="All",
                                                   duration=c(3,100),
                                                   failRate=log(2)/c(9,18),
                                                   hr=c(.9,.6),
                                                   dropoutRate=rep(.001,2)),
                          ratio=1,                # Experimental:Control randomization ratio
                          events = NULL, # Events at analyses
                          analysisTimes = NULL,   # Times of analyses
                          rho,
                          gamma,
                          tau =  rep(-1, length(rho)),
                          approx = "asymptotic"
){

  weight <- get_combo_weight(rho, gamma, tau)

  info <- lapply(weight, function(x){
    x <- eval(parse(text = x))
    gs_info_wlr(enrollRates, failRates, ratio, events = events, analysisTimes = analysisTimes, weight = x)
  })

  info <- dplyr::bind_rows(info, .id = "test")
  info$test <- as.numeric(info$test)
  info
}
