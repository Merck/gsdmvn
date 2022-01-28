gs_design_combo2 <- function(
  enrollRates,
  failRates,
  fh_test,
  ratio = 1,
  alpha = 0.025,
  beta = 0.2,
  binding = FALSE,
  upper = gs_b,
  upar = c(3,2,1),
  lower = gs_b,
  lpar = c(-1, 0, 1),
  algorithm = mvtnorm::GenzBretz(maxpts = 1e5, abseps = 1e-5),
  n_upper_bound = 1e3,
  ...){
  
  # --------------------------------------------- #
  #     check input values                        #
  # --------------------------------------------- #
  # Currently only support user defined lower and upper bound
  stopifnot( identical(upper, gs_b) | identical(upper, gs_spending_combo) )
  stopifnot( identical(lower, gs_b) | identical(lower, gs_spending_combo) )
  
  # --------------------------------------------- #
  #     get the number of analysis                #
  # --------------------------------------------- #
  n_analysis <- length(unique(fh_test$Analysis))
  
  # --------------------------------------------- #
  #     obtain utilities                          #
  # --------------------------------------------- #
  utility <- gs_utility_combo(enrollRates = enrollRates,
                              failRates = failRates,
                              fh_test = fh_test,
                              ratio = ratio,
                              algorithm = algorithm, 
                              ...)
  
  info     <- utility$info_all
  info_fh  <- utility$info
  theta_fh <- utility$theta
  corr_fh  <- utility$corr
  
  # Information Fraction
  if(n_analysis == 1){
    min_info_frac <- 1
  }else{
    info_frac <- tapply(info$info0, info$test, function(x) x / max(x))
    min_info_frac <- apply(do.call(rbind, info_frac), 2, min)
  }
  
  
  # Function to calculate power
  foo <- function(n, beta, ...){
    
    # Probability Cross Boundary
    prob <- gs_prob_combo(upper_bound = bound$upper,
                          lower_bound = bound$lower,
                          analysis = info_fh$Analysis,
                          theta = theta_fh * sqrt(n),
                          corr = corr_fh,
                          algorithm = algorithm, ...)
    
    max(subset(prob, Bound == "Upper")$Probability) - (1 - beta)
  }
  
  # Find sample isze and bound
  n <- max(info$N)
  n0 <- 0
  while( (abs(n - n0)) > 1e-2){
    # print(n)
    n0 <- n
    
    # Obtain spending function
    bound <- gs_bound(alpha = upper(upar, min_info_frac),
                      beta = lower(lpar, min_info_frac),
                      analysis = info_fh$Analysis,
                      theta = theta_fh * sqrt(n),
                      corr = corr_fh,
                      binding_lower_bound = binding,
                      algorithm = algorithm,
                      alpha_bound = identical(upper, gs_b),
                      beta_bound = identical(lower, gs_b),
                      ...)
    
    
    n <- uniroot(foo, c(1, n_upper_bound), extendInt = "yes", beta = beta, ...)$root
    
  }
  
  
  # Probability Cross Boundary
  prob <- gs_prob_combo(upper_bound = bound$upper,
                        lower_bound = bound$lower,
                        analysis = info_fh$Analysis,
                        theta = theta_fh * sqrt(n),
                        corr = corr_fh,
                        algorithm = algorithm, ...)
  
  # Probability Cross Boundary under Null
  prob_null <- gs_prob_combo(upper_bound = bound$upper,
                             lower_bound = if(binding){bound$lower}else{rep(-Inf, nrow(bound))},
                             analysis = info_fh$Analysis,
                             theta = rep(0, nrow(info_fh)),
                             corr = corr_fh,
                             algorithm = algorithm, ...)
  
  if(binding == FALSE){
    prob_null$Probability[prob_null$Bound == "Lower"] <- NA
  }
  
  prob$Probability_Null <- prob_null$Probability
  
  # Prepare output
  db <- merge(
    data.frame(Analysis = 1:(nrow(prob)/2), prob, Z = unlist(bound)),
    #tibble::tibble(Analysis = 1:(nrow(prob)/2), prob, Z = unlist(bound)),
    info_fh %>% # unique(info_fh[, c("Analysis", "Time", "N", "Events")])
      tibble::as_tibble() %>% 
      select(Analysis, Time, N, Events) %>% 
      unique()
  )
  
  # update sample size and events
  db$Events <- db$Events * n / max(db$N)
  db$N <- db$N * n / max(db$N)
  
  # db[order(db$Bound, decreasing = TRUE), c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability", "Probability_Null")]
  out <- db[order(db$Bound, decreasing = TRUE), c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability", "Probability_Null")]
  
  out_H1 <- out[, c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability")]
  out_H1$hypothesis <- rep("H1", n_analysis)
  out_H1$`Nominal p` <- pnorm(out_H1$Z * (-1))
  
  out_H0 <- out[, c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability_Null")]
  out_H0$Probability <- out_H0$Probability_Null
  out_H0 <- out_H0[, c("Analysis", "Bound", "Time", "N", "Events", "Z", "Probability")]
  out_H0$hypothesis <- rep("H0", n_analysis)
  out_H0$`Nominal p` <- pnorm(out_H0$Z * (-1))
  
  out <- tibble::tibble(Analysis = c(out_H1$Analysis, out_H0$Analysis),
                        Bound = c(out_H1$Bound, out_H0$Bound),
                        Time = c(out_H1$Time, out_H0$Time),
                        N = c(out_H1$N, out_H0$N),
                        Events = c(out_H1$Events, out_H0$Events),
                        Z = c(out_H1$Z, out_H0$Z),
                        Probability = c(out_H1$Probability, out_H0$Probability),
                        hypothesis = c(out_H1$hypothesis, out_H0$hypothesis),
                        `Nominal p` = c(out_H1$`Nominal p`, out_H0$`Nominal p`))
  #out <- union(out_H1, out_H0)
  
  #out$hypothesis <- rep(c("H1", "H0"), each = 2*n_analysis)
  #out$`Nominal p` <- pnorm(out$Z * (-1))
  
  # out <- union(out %>% select(-Probability_Null),
  #             out %>% select(-Probability) %>% rename(Probability = Probability_Null))
  # out <- out %>%
  #   mutate(hypothesis = rep(c("H1", "H0"), each = 2*n_analysis),
  #          "Nominal p" = pnorm(Z * (-1)))
  
  #out$hypothesis = rep(c("H1", "H0"), each = 2*n_analysis)
  #out[["Nominal p"]] = pnorm(out$Z * (-1))
  
  return(out)
}
