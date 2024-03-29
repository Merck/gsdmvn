% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_info_wlr.R
\name{gs_info_wlr}
\alias{gs_info_wlr}
\title{Information and effect size for Weighted Log-rank test}
\usage{
gs_info_wlr(
  enrollRates = tibble::tibble(Stratum = "All", duration = c(2, 2, 10), rate = c(3, 6,
    9)),
  failRates = tibble::tibble(Stratum = "All", duration = c(3, 100), failRate =
    log(2)/c(9, 18), hr = c(0.9, 0.6), dropoutRate = rep(0.001, 2)),
  ratio = 1,
  events = NULL,
  analysisTimes = NULL,
  weight = wlr_weight_fh,
  approx = "asymptotic"
)
}
\arguments{
\item{enrollRates}{enrollment rates}

\item{failRates}{failure and dropout rates}

\item{ratio}{Experimental:Control randomization ratio}

\item{events}{Targeted minimum events at each analysis}

\item{analysisTimes}{Targeted minimum study duration at each analysis}

\item{weight}{weight of weighted log rank test
\itemize{
\item \code{"1"}=unweighted,
\item \code{"n"}=Gehan-Breslow,
\item \code{"sqrtN"}=Tarone-Ware,
\item \code{"FH_p[a]_q[b]"}= Fleming-Harrington with p=a and q=b
}}

\item{approx}{approximate estimation method for Z statistics
\itemize{
\item \code{"event driven"} = only work under proportional hazard model with log rank test
\item \code{"asymptotic"}
}}
}
\value{
a \code{tibble} with columns \code{Analysis, Time, N, Events, AHR, delta, sigma2, theta, info, info0.}
\code{info, info0} contains statistical information under H1, H0, respectively.
For analysis \code{k}, \code{Time[k]} is the maximum of \code{analysisTimes[k]} and the expected time
required to accrue the targeted \code{events[k]}.
\code{AHR} is expected average hazard ratio at each analysis.
}
\description{
Based on piecewise enrollment rate, failure rate, and dropout rates computes
approximate information and effect size using an average hazard ratio model.
}
\details{
The \code{AHR()} function computes statistical information at targeted event times.
The \code{tEvents()} function is used to get events and average HR at targeted \code{analysisTimes}.
}
