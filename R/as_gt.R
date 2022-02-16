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

#' This is the function to format the bounds summary table into gt style.
#'
#' @param x an object returned by \code{summary_bound}
#' @param title a string to specify the title of the gt table
#' @param subtitle a string to specify the subtitle of the gt table
#' @param colname_spanner a string to specify the spanner of the gt table
#' @param colname_spannersub a vector of strings to specify the spanner details of the gt table
#' @param footnote a list containing \code{content}, \code{location}, and \code{attr}. 
#'                 the \code{content} is a vector of string to specify the footnote text;
#'                 the \code{location} is a vector of string to specify the locations to put the superscript of the footnote index;
#'                 the \code{attr} is a vector of string to specify the attributes of the footnotes, e.g., c("colname", "title", "subtitle", "analysis", "spanner");
#'                 users can use the functions in the \code{gt} package to custom themselves.
#'
#' @return a gt table summarizing the bounds table in group sequential designs
#' @export
#'
#' @examples 
#' x_combo <- gsdmvn::gs_design_combo2(
#'   ratio = 1, 
#'   alpha = 0.025, 
#'   beta = 0.2, 
#'   enrollRates = tibble::tibble(Stratum = "All", duration = 12, rate = 500/12),
#'   failRates = tibble::tibble(Stratum = "All", duration = c(4, 100), 
#'                              failRate = log(2) / 15, hr = c(1, .6), dropoutRate = .001), 
#'   fh_test = fh_test,
#'   upper = gsdmvn::gs_spending_combo,
#'   upar = list(sf = gsDesign::sfLDOF, total_spend = 0.025),
#'   lower = gsdmvn::gs_spending_combo,
#'   lpar = list(sf = gsDesign::sfLDOF, total_spend = 0.2))
#'   summary_bound(x_combo2, method = "COMBO") %>% 
#' as_gt(
#'   title = "Summary of the Crossing Probability",
#'   subtitle = "by Using gs_design_wlr",
#'   colname_spanner = "Cumulative boundary crossing probability",
#'   colname_spannersub = c("Alternate hypothesis", "Null hypothesis"),
#'   footnote = list(content = c("gs_design_wlr is a function in gsdmvn.",
#'                               "EF is the events fraction, i.e., Events/N."),
#'                   location = c(NA, NA),
#'                   attr = c("subtitle", "analysis")))
#'                           
as_gt <- function(
  x,
  title = NULL,
  subtitle = NULL,
  colname_spanner = "Cumulative boundary crossing probability",
  colname_spannersub = c("Alternate hypothesis", "Null hypothesis"),
  footnote = NULL
){
  method <- class(x)[class(x) %in% c("ahr", "wlr", "combo")]
  # --------------------------------------------- #
  #     set defaults                              #
  # --------------------------------------------- #
  # set different default title to different methods
  if(method == "ahr" && is.null(title)){
    title <- "Bound summary for gs_design_ahr"
  }
  if(method == "wlr" && is.null(title)){
    title <- "Bound summary for gs_design_wlr"
  }
  if(method == "combo" && is.null(title)){
    title <- "Bound summary for gs_design_combo"
  }
  
  # set different default subtitle to different methods
  if(method == "ahr" && is.null(subtitle)){
    subtitle <- "AHR approximations of ~HR at bound"
  }
  if(method == "wlr" && is.null(subtitle)){
    subtitle <- "WLR approximation of ~wHR at bound"
  }
  if(method == "combo" && is.null(subtitle)){
    subtitle <- "Max Combo approximation"
  }
  
  # set different default footnotes to different methods
  if(method == "ahr" && is.null(footnote)){
    footnote <- list(content = "approximate hazard ratio to cross bound.",
                     location = "~HR at bound",
                     attr = "colname")
  }
  if(method == "wlr" && is.null(footnote)){
    footnote <- list(content = c("approximate weighted hazard ratio to cross bound.", "wAHR is the weighted AHR."),
                     location = c("~wHR at bound", NA),
                     attr = c("colname", "analysis"))
  }
  if(method == "combo" && is.null(footnote)){
    footnote <- list(content = "EF is event fraction. AHR  is under regular weighted log rank test.",
                     location = NA,
                     attr = "analysis")
  }
  
  # add spanner 
  x <- x %>% 
    dplyr::group_by(Analysis) %>%
    gt::gt() %>%
    gt::tab_spanner(
      columns = all_of(colname_spannersub),
      label = colname_spanner) %>%
    gt::fmt_number(
      columns = all_of(seq(3, ifelse(method %in% c("ahr", "wlr"), 6, 5), by = 1)), #3:6
      decimals = 4) %>%
    gt::tab_header(
      title = title,
      subtitle = subtitle)
  
  # add footnotes
  if(!is.null(footnote$content)){
    for (i in 1:length(footnote$content)) {
      # if the footnotes is added on the colnames
      if(footnote$attr[i] == "colname"){
        x <- x %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_column_labels(columns = footnote$location[i]))
      }
      # if the footnotes is added on the title/subtitle
      if(footnote$attr[i] == "title" || footnote$attr[i] == "subtitle"){
        x <- x %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_title(group = footnote$attr[i]))
      }
      # if the footnotes is added on the analysis summary row, which is a grouping variable, i.e., Analysis
      if(footnote$attr[i] == "analysis"){
        x <- x %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_row_groups(groups = dplyr::starts_with("Analysis")))
      }
      # if the footnotes is added on the column spanner
      if(footnote$attr[i] == "spanner"){
        x <- x %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_column_spanners(spanners = colname_spanner)
          )
      }
    }
  }
  return(x)
}