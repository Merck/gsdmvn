gs_design_npe2 <- function(
  theta = .1, theta1 = NULL, 
  info = 1, info0 = NULL, info1 = NULL,
  alpha = 0.025, beta = .1, binding = FALSE,
  upper = gs_b, upar = qnorm(.975), test_upper = TRUE,
  lower = gs_b, lpar= -Inf, test_lower = TRUE,
  r = 18, tol = 1e-6,
  calc_h0_prob = FALSE){
  #######################################################################################
  # WRITE INPUT CHECK TESTS AND RETURN APPROPRIATE ERROR MESSAGES
  # info should be a scalar or vector of positive increasing values
  # info0, info1 should be NULL or of the same form as info
  # theta should be a scalar or vector of real values; if vector, same length as info
  # theta0, theta1 should be NULL or same form and length as theta
  # test_upper and test_lower should be logical scalar or vector; if vector same length as info
  # alpha and beta should be scalars with 0 < alpha < 1 - beta < 1
  
  # --------------------------------------------- #
  #     check info, info0, info1                  #
  # --------------------------------------------- #
  if (!is.vector(info, mode = "numeric")) stop("gs_design_npe(): info must be specified numeric vector")
  K <- length(info)
  if (is.null(info0)) info0 <- info
  if (is.null(info1)) info1 <- info
  if (!is.vector(info0, mode = "numeric")) stop("gs_design_npe(): info0 must be specified numeric vector or NULL")
  if (!is.vector(info1, mode = "numeric")) stop("gs_design_npe(): info1 must be specified numeric vector or NULL")
  if (length(info1) != length(info) || length(info0) != length(info) ) stop("gs_design_npe(): length of info, info0, info1 must be the same")
  if (min(info - lag(info,default = 0) <= 0)) stop("gs_design_npe(): info much be strictly increasing and positive")
  if (min(info0 - lag(info0,default = 0) <= 0)) stop("gs_design_npe(): info0 much be NULL or strictly increasing and positive")
  if (min(info1 - lag(info1,default = 0) <= 0)) stop("gs_design_npe(): info1 much be NULL or strictly increasing and positive")
  
  # --------------------------------------------- #
  #     check theta, theta0, theta1               #
  # --------------------------------------------- #
  if (!is.vector(theta, mode = "numeric")) stop("gs_design_npe(): theta must be a real vector")
  if (length(theta) == 1 && K > 1) theta <- rep(theta, K)
  if (length(theta) != K) stop("gs_design_npe(): if length(theta) > 1, must be same as info")
  if (theta[K] <= 0) stop("gs_design_npe(): final effect size must be > 0")
  if (is.null(theta1)){theta1 <- theta}else if (length(theta1)==1) theta1 <- rep(theta1,K)
  if (!is.vector(theta1, mode = "numeric")) stop("gs_design_npe(): theta1 must be a real vector")
  if (length(theta1) != K) stop("gs_design_npe(): if length(theta1) > 1, must be same as info")
  
  # --------------------------------------------- #
  #     check test_upper & test_lower             #
  # --------------------------------------------- #
  # check the correct spec of test_upper & test_lower
  if (length(test_upper) == 1 && K > 1) test_upper <- rep(test_upper, K)
  if (length(test_lower) == 1 && K > 1) test_lower <- rep(test_lower, K)
  ## Check test_upper and test_lower are logical and correct length
  if (!is.vector(test_upper, mode = "logical") || !is.vector(test_lower, mode = "logical"))
    stop("gs_design_npe(): test_upper and test_lower must be logical")
  if (!(length(test_upper) == 1 || length(test_upper) == K))
    stop("gs_design_npe(): test_upper must be length 1 or same length as info")
  if (!(length(test_lower) == 1 || length(test_lower) == K))
    stop("gs_design_npe(): test_lower must be length 1 or same length as info")
  ## Check that final test_upper value is TRUE
  if (!dplyr::last(test_upper)) stop("gs_design_npe(): last value of test_upper must be TRUE")
  
  # --------------------------------------------- #
  #     check alpha & beta                        #
  # --------------------------------------------- #
  ## Check alpha and beta numeric, scalar, 0 < alpha < 1 - beta
  if (!is.numeric(alpha)) stop("gs_design_npe(): alpha must be numeric")
  if (!is.numeric(beta)) stop("gs_design_npe(): beta must be numeric")
  if (length(alpha) != 1 || length(beta) != 1) stop("gs_design_npe(): alpha and beta must be length 1")
  if (alpha <= 0 || 1 - beta <= alpha || beta <= 0) stop("gs_design_npe(): must have 0 < alpha < 1 - beta < 1")
  
  # --------------------------------------------- #
  #     initialization                            #
  # --------------------------------------------- #
  a <- rep(-Inf, K)          # bounds
  b <- rep(Inf, K)
  hgm1_0 <- NULL             # numerical integration grids
  hgm1_1 <- NULL
  upperProb <- rep(NA, K)    # boundary crossing probabilities
  lowerProb <- rep(NA, K)
  
  # --------------------------------------------- #
  #     fixed design                              #
  # --------------------------------------------- #
  # compute fixed sample size for desired power and Type I error.
  minx <- ((qnorm(alpha) / sqrt(info0[K]) + qnorm(beta) / sqrt(info[K])) / theta[K])^2
  # for a fixed design, this is all you need.
  if (K == 1) return(tibble::tibble(
    Analysis = 1,
    Bound = "Upper",
    Z = qnorm(1 - alpha),
    Probability = 1 - beta,
    theta = theta,
    info = info * minx,
    info0 = info0 * minx)
  )
  
  # find an interval for information inflation to give correct power
  minpwr <- gs_power_npe2(
    theta = theta, theta1 = theta1,
    info = info * minx, info1 = info * minx, info0 = info0 * minx,
    binding = binding,
    upper = upper, upar = upar, test_upper = test_upper,
    lower = lower, lpar = lpar, test_lower = test_lower,
    r = r, tol = tol
  ) %>% 
    filter(hypothesis == "H1" & Bound == "Upper" & Analysis == K) %>% 
    select(Probability) %>% 
    unlist() %>% 
    as.numeric()  # $Probability[K]
  
  # --------------------------------------------- #
  #     FOLLOWING IS PAINFUL                      #
  #       AND SHOULD NEVER BE NEEDED              #
  #     BUT IF IT IS NEEDED,                      #
  #       IT TELLS YOU WHAT WENT WRONG!           #
  #     NEED TO BRACKET TARGETED POWER            #
  #       BEFORE ROOT FINDING                     #
  # --------------------------------------------- #
  ## Ensure minx gives power < 1 - beta and maxx gives power > 1 - beta
  if (minpwr < 1 - beta){
    maxx <- 1.05 * minx
    ## Ensure maxx is sufficient information inflation to overpower
    err <- 1
    for(i in 1:10){
      maxpwr <- gs_power_npe2(
        theta = theta, theta1 = theta1,
        info = info * maxx, info1 = info * maxx, info0 = info0 * maxx,
        binding = binding,
        upper = upper, upar = upar, test_upper = test_upper, 
        lower = lower, lpar= lpar, test_lower = test_lower,
        r = r, tol = tol
      )%>% 
        filter(hypothesis == "H1" & Bound == "Upper" & Analysis == K) %>% 
        select(Probability) %>% 
        unlist() %>% 
        as.numeric() #$Probability[K]
      
      if (1  - beta > maxpwr){
        minx <- maxx
        maxx <- 1.05 * maxx
      }else{
        err <- 0
        break
      }
    }
    if (err) stop("gs_design_npe: could not inflate information to bracket power before root finding")
  }else{
    maxx <- minx
    minx <- maxx / 1.05
    err <- 1
    for(i in 1:10){
      if (1  - beta < gs_power_npe2(
        theta = theta, theta1 = theta1,
        info = info * minx, info1 = info1 * minx, info0 = info0 * minx,
        binding = binding,
        upper=upper, lower=lower, upar = upar, lpar= lpar,
        test_upper = test_upper, test_lower = test_lower,
        r = r, tol = tol
      ) %>% 
      filter(hypothesis == "H1" & Bound == "Upper" & Analysis == K) %>% 
      select(Probability) %>% 
      unlist() %>% 
      as.numeric() #$Probability[K]
      ){
        maxx <- minx
        minx <- minx / 1.05}else{
          err <- 0
          break
        }
    }
    if (err) stop("gs_design_npe: could not deflate information to bracket targeted power before root finding")
  }
  
  # --------------------------------------------- #
  #     EITHER TARGETED POWER NOW BRACKETED       #
  #                   OR                          #
  #     ERROR MESSAGE HAS BEEN RETURNED           #
  #     AND WE CAN ACTUALLY GO ON TO FIND THE ROOT#
  # --------------------------------------------- #
  ## Use root finding with the above function to find needed sample size inflation
  # Now we can solve for the inflation factor for the enrollment rate to achieve the desired power
  res <- try(
    uniroot(errbeta, lower = minx, upper = maxx,
            theta = theta, theta1 = theta1, K = K, beta = beta,
            info = info, info1 = info1, info0 = info0, binding = binding,
            Zupper = upper, Zlower = lower, upar = upar, lpar = lpar,
            test_upper = test_upper, test_lower = test_lower,
            r = r, tol = tol)
  )
  if(inherits(res, "try-error")){stop("gs_design_npe: Sample size solution not found")}
  
  ## Update targeted info, info0 based on inflation factor and return a tibble with
  ## bounds, targeted information, and boundary crossing probabilities at each analysis
  # return(gs_power_npe(theta = theta, theta1 = theta1,
  #                     info = info * res$root, info1 = info1 * res$root, info0 = info0 * res$root,
  #                     binding = binding,
  #                     upper=upper, lower=lower, upar = upar, lpar= lpar,
  #                     test_upper = test_upper, test_lower = test_lower,
  #                     r = r, tol = tol))
  
  if(calc_h0_prob == TRUE){
    return(gs_power_npe2(
      theta = theta, theta1 = theta1,
      info = info * res$root, info1 = info1 * res$root, info0 = info0 * res$root,
      binding = binding,
      upper = upper, upar = upar, test_upper = test_upper,
      lower = lower, lpar = lpar, test_lower = test_lower,
      r = r, tol = tol,
      calc_h0_prob = TRUE))
  }else{
    return(gs_power_npe2(
      theta = theta, theta1 = theta1,
      info = info * res$root, info1 = info1 * res$root, info0 = info0 * res$root,
      binding = binding,
      upper = upper, upar = upar, test_upper = test_upper,
      lower = lower, lpar = lpar, test_lower = test_lower,
      r = r, tol = tol))
  }
  
}
## Create a function that uses gs_power_npe to compute difference from targeted power
## for a given sample size inflation factor
errbeta <- function(x = 1, K = 1, beta = .1, theta = .1, theta1 = .1, info = 1, info1 = 1, info0 = 1, binding = FALSE,
                    Zupper = gs_b, Zlower = gs_b, upar = qnorm(.975), lpar= -Inf,
                    test_upper = TRUE, test_lower = TRUE,
                    r = 18, tol = 1e-6){
  return(1 -  beta -
           gs_power_npe2(theta = theta, theta1 = theta1,
                         info = info * x, info1 = info1 * x, info0 = info0 * x, binding = binding,
                         upper = Zupper, lower = Zlower, upar = upar, lpar= lpar,
                         test_upper = test_upper, test_lower = test_lower,
                         r = r, tol = tol
           )%>% 
           filter(hypothesis == "H1" & Bound == "Upper" & Analysis == K) %>% 
           select(Probability) %>% 
           unlist() %>% 
           as.numeric()) # $Probability[K])
}