summary_bound <- function(
  x,
  method = NULL,
  gt_style = TRUE, # if gt_style = TRUE, the output a gt table
  title = NULL,
  subtitle = NULL,
  colname_spanner = "Cumulative boundary crossing probability",
  colname_spannersub = c("Alternate hypothesis", "Null hypothesis"),
  analysis_vars = NULL,
  analysis_decimals = NULL,
  bound_names = c("Efficacy", "Futility"),
  footnote = NULL
){
  
  # --------------------------------------------- #
  #     set defaults                              #
  # --------------------------------------------- #
  # set different default title to different methods
  if(method == "AHR" && is.null(title)){
    title <- "Bound summary for gs_design_ahr"
  }
  if(method == "WLR" && is.null(title)){
    title <- "Bound summary for gs_design_wlr"
  }
  if(method == "COMBO" && is.null(title)){
    title <- "Bound summary for gs_design_combo"
  }
  
  # set different default subtitle to different methods
  if(method == "AHR" && is.null(subtitle)){
    subtitle <- "AHR approximations of ~HR at bound"
  }
  if(method == "WLR" && is.null(subtitle)){
    subtitle <- "WLR approximation of ~wHR at bound"
  }
  if(method == "COMBO" && is.null(subtitle)){
    subtitle <- "COMBO approximations of XXX at bound"
  }
  
  # set different default footnotes to different methods
  if(method == "AHR" && is.null(footnote)){
    footnote <- list(content = "approximate hazard ratio to cross bound.",
                     location = "~HR at bound",
                     attr = "colname")
  }
  if(method == "WLR" && is.null(footnote)){
    footnote <- list(content = "approximate weighted hazard ratio to cross bound.",
                     location = "~wHR at bound",
                     attr = "colname")
  }
  if(method == "COMBO" && is.null(footnote)){
    footnote <- list(content = "XXX",
                     location = "Analysis",
                     attr = "colname")
  }
  
  # set different default summary header to different methods
  if((method == "AHR" || method == "WLR") && is.null(analysis_vars)){
    analysis_vars <- c("Time", "N", "Events", "AHR", "IF")
    analysis_decimals <- c(1, 1, 1, 2, 2)
  }
  if(method == "COMBO" && is.null(analysis_vars)){
    analysis_vars <- c("Time", "N", "Events")  # TO BE IMPROVED
    analysis_decimals <- c(1, 1, 1)
  }
  
  # calculate the approximate hazard ratio to cross bound
  if(method == "AHR" || method == "WLR"){
    x_bounds <- x$bounds
    K <- max(x$bounds$Analysis)
    temp <- c("Analysis", "Bound", "Nominal p", "~HR at bound", "Alternate hypothesis")
  }
  if(method == "COMBO"){
    x_bounds <- x
    K <- max(x$Analysis)
    temp <- c("Analysis", "Bound", "Nominal p", "Alternate hypothesis")
  }
  
  # class(gs_bounds) <- c("bound", class(x_bounds))
  
  # output <- bound.table(x_bounds,
  #                       analysis_vars = analysis_vars,
  #                       analysis_decimals = analysis_decimals)
  # Make a table with one row per analysis with Analysis, Time, N, Events, AHR, IF
  if(method == "AHR" || method == "WLR"){
    analyses <- x_bounds %>%
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
  
  # Compute null hypothesis table and add information fraction (IF)
  # This information should probably be added already in output bounds from original functions
  # TBD This should be moved to individual functions so it will be expected in input
  # y0 <- gsdmvn::gs_power_npe(
  #   theta = rep(0, nrow(analyses)),
  #   info = analyses$info0,
  #   info0 = analyses$info0,
  #   # Just specify Z-values for bounds
  #   upper  = gsdmvn::gs_b,
  #   lower = gsdmvn::gs_b,
  #   # Note that we must have an upper and a lower bound for each analysis
  #   # (any bound not to be used will be infinite)
  #   upar = (bounds %>% filter(Bound == "Upper"))$Z,
  #   lpar = (bounds %>% filter(Bound == "Lower"))$Z
  # ) %>% filter(abs(Z) < Inf)
  
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
      columns = all_of(colname_spannersub),
      label = colname_spanner) %>%
    gt::fmt_number(
      columns = all_of(tmp), #3:6
      decimals = 4) %>%
    gt::tab_header(
      title = title,
      subtitle = subtitle)
  
  # add footnotes
  if(!is.null(footnote$content)){
    for (i in 1:length(footnote$content)) {
      # if the footnotes is added on the colnames
      if(footnote$attr[i] == "colname"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_column_labels(columns = footnote$location[i]))
      }
      # if the footnotes is added on the title/subtitle
      if(footnote$attr[i] == "title" || footnote$attr[i] == "subtitle"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_title(group = footnote$attr[i]))
      }
      # if the footnotes is added on the analysis summary row, which is a grouping variable, i.e., Analysis
      if(footnote$attr[i] == "analysis"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_group(groups = starts_with("Analysis")))
      }
      # if the footnotes is added on the column spanner
      if(footnote$attr[i] == "spanner"){
        output <- output %>% 
          gt::tab_footnote(
            footnote = footnote$content[i],
            locations = gt::cells_column_spanners(spanners = colname_spanner)
          )
      }
    }
  }
  
  return(output)
}