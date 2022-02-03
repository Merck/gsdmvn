#' @title Group sequential bound summary table
#' @description
#' Produce a concise table summarizing bounds and operating characteristics under the null and an alternative hypothesis for a group sequential design.
#' Initially designed to work with the \pkg{gt} package.
#' Works with designs and boundary summaries produced using
#' \code{\link{gs_design_ahr}}, \code{\link{gs_power_ahr}},
#' \code{\link{gs_design_wlr}}, \code{\link{gs_power_wlr}},
#' \code{\link{gs_design_npe}},  or \code{\link{gs_power_npe}}; any design where null hypothesis boundary crossing probabilities can be produced using the
#' \code{\link{gs_design_npe}} function can also be used.
#' Currently, \code{\link{gs_design_combo}} and \code{\link{gs_power_combo}} are not supported.
#'
#' @details
#' Depends on two internal functions: \code{\link{table_ab}} to merge two tables with a one-to-many relationship for printing with \pkg{gt} and
#' \code{\link{rounddf}} for rounding numbers in a data frame.
#' Current version depends on \code{\link{gs_design_npe}} for null hypothesis boundary crossing probability computations,
#' thus excluding the possibility of producing summary tables for the \code{\link{gs_design_combo}} and \code{\link{gs_power_combo}} bound tables.
#'
#' @param bounds Bounds for group sequential design with class `bound`.
#' @param bound_names Names for bounds; default = c("Efficacy", "Futility").
#' @param bound_vars indicates optional columns to be included in boundary table.
#' @param analysis_vars Name(s) related to accrual-related variables; default = c("Time", "N", "Events", "IF", AHR") to be combined with analysis number
#' in a character string summary with one line per analysis.
#' @param analysis_decimals Corresponding number of decimals to display; default = c(1, 1, 1, des2, 2)
#' @param maxinfo0 Denominator for information fraction under null hypothesis;
#' default = NULL in which case largest value in \code{info0} from the input table is used.
#' @param h0_table Default is \code{NULL} in which case null hypothesis boundary crossing probabilities are produced using \code{\link{gs_power_npe}}.
#' This variable is currently not implemented, but intended for future use with \code{\link{gs_design_combo}} and \code{\link{gs_power_combo}}.
#' @return A data frame grouped by `Analysis`.
#' The variable `Analysis`  produced is a character string combining descriptive information for each analysis; this is to be printed once for each analysis.
#' Other variables always produced with one value for each bound for each analysis are:
#' \item{Bound} to indicate which bound (e.g., Efficacy, Futility) is summarized; this is put at the left-hand side of the table,
#' \item{Alternate hypothesis} to provide cumulative boundary crossing probabilities for crossing a bound under the alternate hypothesis (right-hand side of table), and
#' \item{Null hypothesis} to provide cumulative boundary crossing probabilities for crossing a bound under the null hypothesis (right-hand side of table).
#' The variables indicated in input variable \code{bound_vars} are also included after \code{Bound} and before \code{Alternate hypothesis} in the output table.
#' @export
#'
#' @examples
table.bound <- function(bounds,
                        bound_names = c("Efficacy", "Futility"),
                        bound_vars = c("~HR at bound"),
                        analysis_vars = c("Time", "N", "Events", "AHR", "IF"),
                        analysis_decimals = c(1, 1, 1, 2, 2),
                        maxinfo0 = NULL){
  # Default for maxinfo0 is to get it from bound table
  if (is.null(maxinfo0)) maxinfo0 <- max(bounds$info0)
  # Check if bound vars are in input table
  # TBD
  # Check if analysis_vars are in input table
  # TBD
  # Check if analysis_decimals specified for all analysis_vars
  # TBD
  # Make a table with one row per analysis with Analysis, Time, N, Events, AHR, IF
  analyses <- bounds %>%
    group_by(Analysis) %>%
    filter(row_number() == 1) %>%
    mutate(IF = info0 / maxinfo0) %>%
    select(all_of(c("Analysis", analysis_vars, "info0")))

  # Compute null hypothesis table and add information fraction (IF)
  # This information should probably be added already in output bounds from original functions
  # TBD This should be moved to individual functions so it will be expected in input
  y0 <- gs_power_npe(
    theta = rep(0, nrow(analyses)),
    info = analyses$info0,
    info0 = analyses$info0,
    # Just specify Z-values for bounds
    upper  = gs_b,
    lower = gs_b,
    # Note that we must have an upper and a lower bound for each analysis
    # (any bound not to be used will be infinite)
    upar = (bounds %>% filter(Bound == "Upper"))$Z,
    lpar = (bounds %>% filter(Bound == "Lower"))$Z
  ) %>% filter(abs(Z) < Inf)

  # Merge 3 tables: analyses (one to many), alternate hypothesis table, null hypothesis table
  xy <- full_join(bounds %>%
                    rename("Alternate hypothesis" = Probability) %>%
                    mutate("Nominal p" = pnorm(-Z),
                           Bound = recode(Bound, "Upper" = bound_names[1], "Lower" = bound_names[2])) %>%
                    select(all_of(c("Analysis", "Bound", "Nominal p", bound_vars, "Alternate hypothesis"))),
                  y0 %>%
                    mutate("Null hypothesis" = Probability,
                           Bound = recode(Bound, "Upper" = bound_names[1], "Lower" = bound_names[2])) %>%
                    select(all_of(c("Analysis", "Bound", "Null hypothesis"))),
                  by=c("Analysis", "Bound"))

  # Merge 3 tables: 1 line per analysis, alternate hypothesis table, null hypothesis table
  return(table_ab(analyses %>%
                    select(all_of(c("Analysis", analysis_vars))),
                  xy %>% select(all_of(c("Analysis", "Bound", bound_vars, "Nominal p", "Alternate hypothesis", "Null hypothesis"))),
                  decimals = c(0, analysis_decimals),
                  byvar = "Analysis"))
}
