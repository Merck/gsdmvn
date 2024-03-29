% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_prob_combo.R
\name{gs_prob_combo}
\alias{gs_prob_combo}
\title{MaxCombo Group sequential boundary crossing probabilities}
\usage{
gs_prob_combo(
  upper_bound,
  lower_bound,
  analysis,
  theta,
  corr,
  algorithm = GenzBretz(maxpts = 1e+05, abseps = 1e-05),
  ...
)
}
\arguments{
\item{upper_bound}{a numeric vector of upper bound}

\item{lower_bound}{a numeric vector of lower bound}

\item{analysis}{an integer vector of the interim analysis index}

\item{theta}{a numeric vector of effect size under alternative hypothesis}

\item{corr}{a matrix of correlation matrix}

\item{algorithm}{ an object of class \code{\link[mvtnorm]{GenzBretz}},
                    \code{\link[mvtnorm]{Miwa}} or \code{\link[mvtnorm]{TVPACK}}
                    specifying both the algorithm to be used as well as
                    the associated hyper parameters.}

\item{...}{additional parameters transfer to \code{mvtnorm::pmvnorm}}
}
\description{
MaxCombo Group sequential boundary crossing probabilities
}
