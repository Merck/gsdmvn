summary_bound <- function(
  x,
  method = NULL,
  analysis_vars = NULL,
  analysis_decimals = NULL,
  bound_names = c("Efficacy", "Futility"),
  gt_style = FALSE, # if gt_style = TRUE, the output a gt table
  gt_title = NULL,
  gt_subtitle = NULL,
  gt_colname_spanner = "Cumulative boundary crossing probability",
  gt_colname_spannersub = c("Alternate hypothesis", "Null hypothesis"),
  gt_footnote = NULL
){
  
  # get the 
  # (1) bounds table,
  # (2) analysis summary table, 
  # (3) analysis variables to be displayed on the header
  # (4) decimals to be displayed for the analysis variables in (3)
  if(method == "AHR" || method == "WLR"){
    x_bounds <- x$bounds
    x_analysis <- x$analysis
    K <- max(x_analysis$Analysis)
    analysis_vars <- c("Time", "N", "Events", "AHR", "IF")
    analysis_decimals <- c(1, 1, 1, 2, 2)
    temp <- c("Analysis", "Bound", "Nominal p", "~HR at bound", "Alternate hypothesis")
  }
  if(method == "COMBO"){
    x_bounds <- x
    K <- max(x$Analysis)
    analysis_vars <- c("Time", "N", "Events")  # TO BE IMPROVED
    analysis_decimals <- c(1, 1, 1)
    temp <- c("Analysis", "Bound", "Nominal p", "Alternate hypothesis")
  }
  
  # set the analysis summary header
  if(method == "AHR" || method == "WLR"){
    analyses <- x_analysis %>%
      group_by(Analysis) %>%
      filter(row_number() == 1) %>%
      select(all_of(c("Analysis", analysis_vars, "info0")))
  }
  if(method == "COMBO"){
    analyses <- x_bounds %>%
      group_by(Analysis) %>%
      filter(row_number() == 1) %>%
      select(all_of(c("Analysis", analysis_vars)))
  }
  
  # Merge 3 tables: analyses (one to many), alternate hypothesis table, null hypothesis table
  xy <- full_join(
    # a table under alternative hypothesis
    x_bounds %>% 
      filter(hypothesis == "H1") %>% 
      rename("Alternate hypothesis" = Probability) %>%
      mutate(Bound = recode(Bound, "Upper" = bound_names[1], "Lower" = bound_names[2])
            # change Upper -> bound_names[1], e.g., Efficacy
            # change Lower -> bound_names[2], e.g., Futility
            ) %>%
      select(all_of(temp)), 
      #select(all_of(c("Analysis", "Bound", "Nominal p", "~HR at bound", "Alternate hypothesis"))),
    
    # a table under null hypothesis
    x_bounds %>% 
      filter(hypothesis == "H0") %>%
      mutate("Null hypothesis" = Probability,
             Bound = recode(Bound, "Upper" = bound_names[1], "Lower" = bound_names[2])) %>%
      select(all_of(c("Analysis", "Bound", "Null hypothesis"))),
    
    by = c("Analysis", "Bound"))
  
  # Merge 3 tables: 1 line per analysis, alternate hypothesis table, null hypothesis table
  # if the method is AHR
  if(method == "AHR"){
    # header
    analysis_summary_header <- analyses %>%
      select(all_of(c("Analysis", analysis_vars)))
    # bound details
    bound_summary_detail <- xy %>%
      select(all_of(c("Analysis", "Bound", "~HR at bound", 
                      "Nominal p", "Alternate hypothesis", "Null hypothesis")))
  }
  # if the method is WLR, change AHR to wAHR
  if(method == "WLR"){
    # header
    analysis_summary_header <- analyses %>%
      select(all_of(c("Analysis", analysis_vars)))
    if("AHR" %in% analysis_vars){
      analysis_summary_header <- analysis_summary_header %>% rename(wAHR = AHR)
    }
    # bound details
    bound_summary_detail <- xy %>%
      select(all_of(c("Analysis", "Bound", "~HR at bound", 
                      "Nominal p", "Alternate hypothesis", "Null hypothesis"))) %>% 
      rename("~wHR at bound" = "~HR at bound")
  }
  # if the method is COMBO, remove the column of "~HR at bound", and remove AHR from header 
  if(method == "COMBO"){
    # header
    analysis_summary_header <- analyses %>%
      select(all_of(c("Analysis", analysis_vars)))
    # bound details
    bound_summary_detail <- xy %>%
      select(all_of(c("Analysis", "Bound", 
                      "Nominal p", "Alternate hypothesis", "Null hypothesis"))) 
  }
  
  output <- table_ab(
    # A data frame to be show as the summary header
    # It has only ONE record for each value of `byvar`
    table_a = analysis_summary_header,
    # A data frame to be shown as the listing details
    # It has >= 1 records for each value of `byvar`
    table_b = bound_summary_detail,
    decimals = c(0, analysis_decimals),
    byvar = "Analysis")
  
  if(gt_style == FALSE){
    return(output)
  }
  
  
  
  # --------------------------------------------- #
  #     the following discusses when              #
  #          gt_style = TRUE                      #
  # --------------------------------------------- #
  # --------------------------------------------- #
  #     set defaults                              #
  # --------------------------------------------- #
  # set different default title to different methods
  if(method == "AHR" && is.null(gt_title)){
    gt_title <- "Bound summary for gs_design_ahr"
  }
  if(method == "WLR" && is.null(gt_title)){
    gt_title <- "Bound summary for gs_design_wlr"
  }
  if(method == "COMBO" && is.null(gt_title)){
    gt_title <- "Bound summary for gs_design_combo"
  }
  
  # set different default subtitle to different methods
  if(method == "AHR" && is.null(gt_subtitle)){
    gt_subtitle <- "AHR approximations of ~HR at bound"
  }
  if(method == "WLR" && is.null(gt_subtitle)){
    gt_subtitle <- "WLR approximation of ~wHR at bound"
  }
  if(method == "COMBO" && is.null(gt_subtitle)){
    gt_subtitle <- "COMBO approximations of XXX at bound"
  }
  
  # set different default footnotes to different methods
  if(method == "AHR" && is.null(gt_footnote)){
    gt_footnote <- list(content = "approximate hazard ratio to cross bound.",
                        location = "~HR at bound",
                        attr = "colname")
  }
  if(method == "WLR" && is.null(gt_footnote)){
    gt_footnote <- list(content = "approximate weighted hazard ratio to cross bound.",
                        location = "~wHR at bound",
                        attr = "colname")
  }
  if(method == "COMBO" && is.null(gt_footnote)){
    gt_footnote <- list(content = "XXX",
                        location = "Analysis",
                        attr = "colname")
  }
  
  # add spanner  -- TO BE IMPROVED
  if(method %in% c("AHR", "WLR")){
    tmp <- 3:6
  }else{
    tmp <- 3:5
  }
  output <- output %>% 
    dplyr::group_by(Analysis) %>%
    gt::gt() %>%
    gt::tab_spanner(
      columns = all_of(gt_colname_spannersub),
      label = gt_colname_spanner) %>%
    gt::fmt_number(
      columns = all_of(tmp), #3:6
      decimals = 4) %>%
    gt::tab_header(
      title = gt_title,
      subtitle = gt_subtitle)
  
  # add footnotes
  if(!is.null(gt_footnote$content)){
    for (i in 1:length(gt_footnote$content)) {
      # if the footnotes is added on the colnames
      if(gt_footnote$attr[i] == "colname"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = gt_footnote$content[i],
            locations = gt::cells_column_labels(columns = gt_footnote$location[i]))
      }
      # if the footnotes is added on the title/subtitle
      if(gt_footnote$attr[i] == "title" || gt_footnote$attr[i] == "subtitle"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = gt_footnote$content[i],
            locations = gt::cells_title(group = gt_footnote$attr[i]))
      }
      # if the footnotes is added on the analysis summary row, which is a grouping variable, i.e., Analysis
      if(gt_footnote$attr[i] == "analysis"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = gt_footnote$content[i],
            locations = gt::cells_group(groups = starts_with("Analysis")))
      }
      # if the footnotes is added on the column spanner
      if(gt_footnote$attr[i] == "spanner"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = gt_footnote$content[i],
            locations = gt::cells_column_spanners(spanners = gt_colname_spanner)
          )
      }
    }
  }
  
  return(output)
}