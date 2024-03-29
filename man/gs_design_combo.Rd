% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_design_combo.R
\name{gs_design_combo}
\alias{gs_design_combo}
\title{Group sequential design using MaxCombo test under non-proportional hazards}
\usage{
gs_design_combo(
  enrollRates,
  failRates,
  fh_test,
  ratio = 1,
  alpha = 0.025,
  beta = 0.2,
  binding = FALSE,
  upper = gs_b,
  upar = c(3, 2, 1),
  lower = gs_b,
  lpar = c(-1, 0, 1),
  algorithm = GenzBretz(maxpts = 1e+05, abseps = 1e-05),
  n_upper_bound = 1000,
  ...
)
}
\arguments{
\item{enrollRates}{enrollment rates}

\item{failRates}{failure and dropout rates}

\item{fh_test}{a data frame to summarize the test in each analysis.
Refer examples for its data structure.}

\item{ratio}{Experimental:Control randomization ratio (not yet implemented)}

\item{alpha}{One-sided Type I error}

\item{beta}{Type II error}

\item{binding}{indicator of whether futility bound is binding; default of FALSE is recommended}

\item{upper}{Function to compute upper bound}

\item{upar}{Parameter passed to \code{upper()}}

\item{lower}{Function to compute lower bound}

\item{lpar}{Parameter passed to \code{lower()}}

\item{algorithm}{ an object of class \code{\link[mvtnorm]{GenzBretz}},
                    \code{\link[mvtnorm]{Miwa}} or \code{\link[mvtnorm]{TVPACK}}
                    specifying both the algorithm to be used as well as
                    the associated hyper parameters.}

\item{n_upper_bound}{a numeric value of upper limit of sample size}

\item{...}{additional parameters transfer to \code{mvtnorm::pmvnorm}}
}
\description{
Group sequential design using MaxCombo test under non-proportional hazards
}
\examples{
\dontrun{

# The example is slow to run

library(dplyr)
library(mvtnorm)
library(gsDesign)

enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 500/12)

failRates <- tibble::tibble(Stratum = "All",
                            duration = c(4, 100),
                            failRate = log(2) / 15,  # median survival 15 month
                            hr = c(1, .6),
                            dropoutRate = 0.001)

fh_test <- rbind( data.frame(rho = 0, gamma = 0, tau = -1,
                             test = 1,
                             Analysis = 1:3,
                             analysisTimes = c(12, 24, 36)),
                  data.frame(rho = c(0, 0.5), gamma = 0.5, tau = -1,
                             test = 2:3,
                             Analysis = 3, analysisTimes = 36)
)

x <- gsDesign::gsSurv( k = 3 , test.type = 4 , alpha = 0.025 ,
                       beta = 0.2 , astar = 0 , timing = c( 1 ) ,
                       sfu = sfLDOF , sfupar = c( 0 ) , sfl = sfLDOF ,
                       sflpar = c( 0 ) , lambdaC = c( 0.1 ) ,
                       hr = 0.6 , hr0 = 1 , eta = 0.01 ,
                       gamma = c( 10 ) ,
                       R = c( 12 ) , S = NULL ,
                       T = 36 , minfup = 24 , ratio = 1 )

# User defined boundary
gs_design_combo(enrollRates,
                failRates,
                fh_test,
                alpha = 0.025,
                beta = 0.2,
                ratio = 1,
                binding = FALSE,                 # test.type = 4 non-binding futility bound
                upar = x$upper$bound,
                lpar = x$lower$bound)

# Boundary derived by spending function
gs_design_combo(enrollRates,
                failRates,
                fh_test,
                alpha = 0.025,
                beta = 0.2,
                ratio = 1,
                binding = FALSE,                 # test.type = 4 non-binding futility bound
                upper = gs_spending_combo,
                upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),   # alpha spending
                lower = gs_spending_combo,
                lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.2),     # beta spending
)
}


}
