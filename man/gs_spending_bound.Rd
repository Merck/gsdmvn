% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs_spending_bound.r
\name{gs_spending_bound}
\alias{gs_spending_bound}
\title{Derive spending bound for group sequential boundary}
\usage{
gs_spending_bound(
  k = 1,
  par = list(sf = gsDesign::sfLDOF, total_spend = 0.025, param = NULL, timing = NULL,
    max_info = NULL),
  hgm1 = NULL,
  theta = 0.1,
  info = 1:3,
  efficacy = TRUE,
  test_bound = TRUE,
  r = 18,
  tol = 1e-06
)
}
\arguments{
\item{k}{analysis for which bound is to be computed}

\item{par}{a list with the following items:
\code{sf} (class spending function),
\code{total_spend} (total spend),
\code{param} (any parameters needed by the spending function \code{sf()}),
\code{timing} (a vector containing values at which spending function is to be evaluated or NULL if information-based spending is used),
\code{max_info} (when \code{timing} is NULL, this can be input as positive number to be used with \code{info} for information fraction at each analysis)}

\item{hgm1}{subdensity grid from h1 (k=2) or hupdate (k>2) for analysis k-1; if k=1, this is not used and may be NULL}

\item{theta}{natural parameter used for lower bound only spending;
represents average drift at each time of analysis at least up to analysis k;
upper bound spending is always set under null hypothesis (theta = 0)}

\item{info}{statistical information at all analyses, at least up to analysis k}

\item{efficacy}{TRUE (default) for efficacy bound, FALSE otherwise}

\item{test_bound}{a logical vector of the same length as \code{info} should indicate which analyses will have a bound}

\item{r}{Integer, at least 2; default of 18 recommended by Jennison and Turnbull}

\item{tol}{Tolerance parameter for convergence (on Z-scale)}
}
\value{
returns a numeric bound (possibly infinite) or, upon failure, generates an error message.
}
\description{
Computes one bound at a time based on spending under given distributional assumptions.
While user specifies \code{gs_spending_bound()} for use with other functions,
it is not intended for use on its own.
Most important user specifications are made through a list provided to functions using \code{gs_spending_bound()}.
Function uses numerical integration and Newton-Raphson iteration to derive an individual bound for a group sequential
design that satisfies a targeted boundary crossing probability.
Algorithm is a simple extension of that in Chapter 19 of Jennison and Turnbull (2000).
}
\section{Specification}{

\if{latex}{
 \itemize{
   \item Set the spending time at analysis.
   \item Compute the cumulative spending at analysis.
   \item Compute the incremental spend at each analysis.
   \item Set test_bound a vector of length k > 1 if input as a single value.
   \item Compute spending for current bound.
   \item Iterate to convergence as in gsbound.c from gsDesign.
   \item Compute subdensity for final analysis in rejection region.
   \item Validate the output and return an error message in case of failure.
   \item Return a numeric bound (possibly infinite).
  }
}
\if{html}{The contents of this section are shown in PDF user manual only.}
}

\references{
Jennison C and Turnbull BW (2000), \emph{Group Sequential
Methods with Applications to Clinical Trials}. Boca Raton: Chapman and Hall.
}
\author{
Keaven Anderson \email{keaven_anderson@merck.com}
}
