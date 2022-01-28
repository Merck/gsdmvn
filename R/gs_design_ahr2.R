gs_design_ahr2 <- function(
  enrollRates = tibble::tibble(Stratum = "All",
                               duration = c(2, 2, 10),
                               rate = c(3, 6, 9)),
  failRates = tibble::tibble(Stratum = "All",
                             duration = c(3, 100),
                             failRate = log(2)/c(9, 18),
                             hr = c(.9, .6),
                             dropoutRate=rep(.001, 2)),
  ratio = 1,             # Experimental:Control randomization ratio
  alpha = 0.025,         # One-sided Type I error
  beta = 0.1,            # NULL if enrollment is not adapted
  IF = NULL,             # relative information fraction timing (vector, if not NULL; increasing to 1)
  analysisTimes = 36,    # Targeted times of analysis or just planned study duration
  binding = FALSE,
  upper = gs_b,
  # Default is Lan-DeMets approximation of
  upar = gsDesign(k = 3, test.type = 1, n.I=c(.25, .75, 1), sfu = sfLDOF, sfupar = NULL)$upper$bound,
  lower = gs_b,
  lpar = c(qnorm(.1), -Inf, -Inf), # Futility only at IA1
  h1_spending = TRUE,
  test_upper = TRUE,
  test_lower = TRUE,
  r = 18,
  tol = 1e-6
){
  # --------------------------------------------- #
  #     check input values                        #
  # --------------------------------------------- #
  msg <- "analysisTimes must be a positive number or positive increasing sequence"
  if (!is.vector(analysisTimes,mode = "numeric")) stop(msg)
  if (min(analysisTimes - dplyr::lag(analysisTimes, def=0))<=0) stop(msg)
  msg <- "gs_design_ahr(): IF must be a positive number or positive increasing sequence on (0, 1] with final value of 1"
  if (is.null(IF)){IF <- 1}
  if (!is.vector(IF,mode = "numeric")) stop(msg)
  if (min(IF - dplyr::lag(IF, def=0))<=0) stop(msg)
  if (max(IF) != 1) stop(msg)
  msg <- "gs_design_ahr() IF and analysisTimes must have the same length if both have length > 1"
  if ((length(analysisTimes)>1) & (length(IF) > 1) & (length(IF) != length(analysisTimes))) stop(msg)
  
  # --------------------------------------------- #
  #     get information at input analysisTimes    #
  # --------------------------------------------- #
  y <- gs_info_ahr(
    enrollRates, 
    failRates,
    ratio = ratio, 
    events = NULL,
    analysisTimes = analysisTimes)
  finalEvents <- y$Events[nrow(y)]
  IFalt <- y$Events / finalEvents
  
  # --------------------------------------------- #
  #     check if IF needed for (any) IA timing    #
  # --------------------------------------------- #
  K <- max(length(analysisTimes), length(IF))
  nextTime <- max(analysisTimes)
  if(length(IF) == 1){
    IF <- IFalt
  }else{
    IFindx <- IF[1:(K-1)]
    for(i in seq_along(IFindx)){
      if(length(IFalt) == 1){
        y <- rbind(
          gsDesign2::tEvents(
            enrollRates, failRates, targetEvents = IF[K - i] * finalEvents, 
            ratio = ratio, interval = c(.01, nextTime)) %>% 
            mutate(theta = -log(AHR), Analysis = K - i),
          y)
      }else if(IF[K-i] > IFalt[K-i]) y[K - i,] <-
          gsDesign2::tEvents(
            enrollRates, failRates, targetEvents = IF[K - i] * finalEvents, ratio = ratio,
            interval = c(.01, nextTime)) %>%
          dplyr::transmute(Analysis = K - i, Time, Events, AHR, theta=-log(AHR), info, info0)
      nextTime <- y$Time[K - i]
    }
  }
  y$Analysis <- 1:K
  y$N <- gsDesign2::eAccrual(x = y$Time, enrollRates = enrollRates)
  if(h1_spending){
    theta1 <- y$theta
    info1 <- y$info
  }else{
    theta1 <- 0
    info1 <- y$info0
  }
  
  # --------------------------------------------- #
  #     get sample size and bounds                #
  #         using gs_design_npe                   #
  # --------------------------------------------- #
  bounds <- gs_design_npe2(
    theta = y$theta, theta1 = theta1,
    info = y$info, info0 = y$info0, info1 = info1,
    alpha = alpha, beta = beta, binding = binding,
    upper = upper, upar = upar, test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol, 
    calc_h0_prob = TRUE) %>%
    # Add Time, Events, AHR, N from gs_info_ahr call above
    full_join(y %>% select(-c(info, info0, theta)), by = "Analysis") %>%
    select(c("Analysis", "Bound", "Time",
             "N", "Events", 
             "Z", "Probability",
             "AHR", "theta", 
             "info", "info0", "IF", "hypothesis")) %>% # "AHR", "theta", "info", "info0")) %>%
    arrange(desc(hypothesis), desc(Bound), Analysis)   # arrange(desc(Bound), Analysis, desc(hypothesis))
  
  bounds$Events <- bounds$Events * bounds$info[K] / y$info[K]
  bounds$N <- bounds$N * bounds$info[K] / y$info[K]
  bounds <- bounds %>% 
    mutate("~HR at bound" = exp(-Z / sqrt(info0)),
           "HR generic (H0)" = exp(-Z / sqrt(info0)),
           "HR generic (H1)" = exp(-Z / sqrt(info)),
           "Nominal p" = pnorm(-Z)) 
  
  return(
    list(enrollRates = enrollRates %>% mutate(rate = rate * bounds$info[K] / y$info[K]),
         failRates = failRates,
         bounds = bounds)
  )
  
}