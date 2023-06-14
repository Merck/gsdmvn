# gsdmvn

## Status [![Lifecycle: superseded](https://img.shields.io/badge/lifecycle-superseded-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#superseded)

gsdmvn is superseded: the functionality was included directly in the gsDesign2 package.
We recommend using [gsDesign2](https://github.com/Merck/gsDesign2) instead.

## Introduction

The goal of **gsdmvn** is to enable fixed or group sequential design
under non-proportional hazards. Piecewise constant enrollment, failure
rates and dropout rates for a stratified population are available to
enable highly flexible enrollment, time-to-event and time-to-dropout
assumptions. Substantial flexibility on top of what is in the gsDesign
package is intended for selecting boundaries. While this work is in
progress, substantial capabilities have been enabled. Comments on
usability and features are encouraged as this is a development version
of the package.

The goal of gsdmvn is to enable group sequential trial design for
time-to-event endpoints under non-proportional hazards assumptions. The
package is still maturing; as the package functions become more stable,
they will likely be included in the
[gsDesign2](https://github.com/Merck/gsDesign2) package.

## Branch specification

-   The `development` branch includes all work under development.
-   The `table_bound` branch is branched from the `development` branch,
    and it targets to get the outputs of `gd_design_ahr()`,
    `gs_power_ahr()`,`gs_desgin_wlr()`, etc., into a well-organized
    form.
-   The `update_futility_bound` branch is branched from the
    `development` branch, and it targets to develop code so one can
    update the futility bound.

## Installation

You can install `gsdmvn` with:

``` r
remotes::install_github("Merck/gsdmvn")
```

## Specifying enrollment and failure rates

This is a basic example which shows you how to solve a common problem.
We assume there is a 4 month delay in treatment effect. Specifically, we
assume a hazard ratio of 1 for 4 months and 0.6 thereafter. For this
example we assume an exponential failure rate and low exponential
dropout rate. The `enrollRates` specification indicates an expected
enrollment duration of 12 months with exponential inter-arrival times.

``` r
library(gsdmvn)
library(gsDesign)
library(gsDesign2)
library(dplyr)
library(knitr)

## basic example code

## Constant enrollment over 12 months
## rate will be adjusted later by gsDesignNPH to get sample size
enrollRates <- tibble::tibble(Stratum = "All", duration = 12, rate = 1)

## 12 month median exponential failure rate in control
## 4 month delay in effect with HR=0.6 after
## Low exponential dropout rate
medianSurv <- 12
failRates <- tibble::tibble(
  Stratum = "All",
  duration = c(4, Inf),
  failRate = log(2) / medianSurv,
  hr = c(1, .6),
  dropoutRate = .001
)
```

The resulting failure rate specification is the following table. As many
rows and strata as needed can be specified to approximate whatever
patterns you wish.

``` r
failRates %>% kable()
```

| Stratum | duration |  failRate |  hr | dropoutRate |
|:--------|---------:|----------:|----:|------------:|
| All     |        4 | 0.0577623 | 1.0 |       0.001 |
| All     |      Inf | 0.0577623 | 0.6 |       0.001 |

Computing a fixed sample size design with 2.5% one-sided Type I error
and 90% power. We specify a trial duration of 36 months with
`analysisTimes`. Since there is a single analysis, we specify an upper
p-value bound of 0.025 with `upar = qnorm(0.975)`. There is no lower
bound which is specified with `lpar = -Inf`.

``` r
design <-
  gs_design_ahr(enrollRates, failRates, upar = qnorm(.975), lpar = -Inf, IF = 1, analysisTimes = 36)
```

The input enrollment rates are scaled to achieve power:

``` r
design$enrollRates %>% kable()
```

| Stratum | duration |     rate |
|:--------|---------:|---------:|
| All     |       12 | 35.05288 |

The failure and dropout rates remain unchanged from what was input:

``` r
design$failRates %>% kable()
```

| Stratum | duration |  failRate |  hr | dropoutRate |
|:--------|---------:|----------:|----:|------------:|
| All     |        4 | 0.0577623 | 1.0 |       0.001 |
| All     |      Inf | 0.0577623 | 0.6 |       0.001 |

Finally, the expected analysis time is in `Time`, sample size `N`,
events required `Events` and bound `Z` are in `design$bounds`. Note that
`AHR` is the average hazard ratio used to calculate the targeted event
counts. The natural parameter (`log(AHR)`) is in theta and corresponding
statistical information under the alternate hypothesis are in `info` and
under the null hypothesis in `info0`.

``` r
design$bounds %>% kable()
```

| Analysis | Bound | Time |        N |   Events |        Z | Probability |       AHR |     theta |     info |    info0 |
|---------:|:------|-----:|---------:|---------:|---------:|------------:|----------:|----------:|---------:|---------:|
|        1 | Upper |   36 | 420.6346 | 311.0028 | 1.959964 |         0.9 | 0.6917244 | 0.3685676 | 76.74383 | 77.75069 |
