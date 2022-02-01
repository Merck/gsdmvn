gs_power_npe2 <- function(
  theta = .1, theta1 = NULL,
  info = 1, info1 = NULL, info0 = NULL,
  binding = FALSE,
  upper = gs_b, upar = qnorm(.975), test_upper = TRUE,
  lower = gs_b, lpar = -Inf, test_lower = TRUE,
  r = 18, tol = 1e-6){
  # --------------------------------------------- #
  #     check & set up parameters                 #
  # --------------------------------------------- #
  K <- length(info)
  if (is.null(info0)) info0 <- info
  if (is.null(info1)) info1 <- info
  if (length(info1) != length(info) || length(info0) != length(info)) stop("gs_power_npe: length of info, info0, info1 must be the same")
  if (length(theta) == 1 && K > 1) theta <- rep(theta, K)
  if (is.null(theta1)){theta1 <- theta}else if (length(theta1)==1) theta1 <- rep(theta1, K)
  if (length(test_upper) == 1 && K > 1) test_upper <- rep(test_upper, K)
  if (length(test_lower) == 1 && K > 1) test_lower <- rep(test_lower, K)
  
  # --------------------------------------------- #
  #     initialization                            #
  # --------------------------------------------- #
  a <- rep(-Inf, K)
  b <- rep(Inf, K)
  hgm1_0 <- NULL
  hgm1_1 <- NULL
  hgm1 <- NULL
  upperProb <- rep(NA, K)
  lowerProb <- rep(NA, K)
  
  # --------------------------------------------- #
  #     calculate crossing probability            #
  #    under the alternative hypothesis           #
  #            i.e., theta != 0                   #
  # --------------------------------------------- #
  for(k in 1:K){
    # compute/update lower bound
    a[k] <- lower(k = k, par = lpar, hgm1 = hgm1_1, info = info1,
                  r = r, tol = tol, test_bound = test_lower,
                  theta = theta1, efficacy = FALSE)
    # compute/update upper bound
    b[k] <- upper(k = k, par = upar, hgm1 = hgm1_0, info = info0,
                  r = r, tol = tol, test_bound = test_upper)
    
    if(k == 1){
      # compute the probability to cross upper/lower bound
      upperProb[1] <- if(b[1] < Inf) {pnorm(b[1], mean = sqrt(info[1]) * theta[1], lower.tail = FALSE)}else{0}
      lowerProb[1] <- if(a[1] > -Inf){pnorm(a[1], mean = sqrt(info[1]) * theta[1])}else{0}
      # update the grids
      hgm1_0 <- h1(r = r, theta = 0,         I = info0[1], a = if(binding){a[1]}else{-Inf}, b = b[1])
      hgm1_1 <- h1(r = r, theta = theta1[1], I = info1[1], a = a[1], b = b[1])
      hgm1   <- h1(r = r, theta = theta[1],  I = info[1],  a = a[1], b = b[1])
    }else{
      # compute the probability to cross upper bound
      upperProb[k] <- if(b[k]< Inf){
        hupdate(r = r, theta = theta[k], I = info[k], a = b[k], b = Inf,
                thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
      }else{0}
      # compute the probability to cross lower bound
      lowerProb[k] <- if(a[k] > -Inf){
        hupdate(r = r, theta = theta[k], I = info[k], a = -Inf, b = a[k],
                thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
      }else{0}
      
      # update the grids
      if(k < K){
        hgm1_0 <- hupdate(r = r, theta = 0, thetam1 = 0, I = info0[k], Im1 = info0[k-1],
                          a = if(binding){a[k]}else{-Inf}, b = b[k], gm1 = hgm1_0)
        hgm1_1 <- hupdate(r = r, theta = theta1[k], thetam1 = theta1[k-1],I = info1[k], Im1 = info1[k-1],
                          a = a[k], b = b[k], gm1 = hgm1_1)
        hgm1   <- hupdate(r = r, theta = theta[k], thetam1 = theta[k-1], I = info[k],  Im1 = info[k-1],
                          a = a[k], b = b[k], gm1 = hgm1)
      }
    }
  }
  
  table_H1 <- tibble::tibble(
    Analysis = rep(1:K, 2),
    Bound = c(rep("Upper", K), rep("Lower", K)),
    Z = c(b, a),
    Probability = c(cumsum(upperProb),
                    cumsum(lowerProb)),
    theta = rep(theta, 2),
    theta1 = rep(theta1, 2),
    IF = rep(info0 / max(info0), 2),
    info = rep(info, 2),
    info0 = rep(info0, 2),
    info1 = rep(info1, 2),
    hypothesis = rep("H1", 2*K))
  
  # --------------------------------------------- #
  #     calculate crossing probability            #
  #       under the null hypothesis               #
  #            i.e., theta == 0                   #
  # --------------------------------------------- #
  # COMPUTE BOUNDS
  # y0 <- gsdmvn::gs_power_npe(
  # theta = rep(0, nrow(analyses)),
  # info = analyses$info0,
  # info0 = analyses$info0,
  # Just specify Z-values for bounds
  # upper  = gsdmvn::gs_b,
  # lower = gsdmvn::gs_b,
  # Note that we must have an upper and a lower bound for each analysis
  # (any bound not to be used will be infinite)
  # upar = (bounds %>% filter(Bound == "Upper"))$Z,
  # lpar = (bounds %>% filter(Bound == "Lower"))$Z
  # ) %>% filter(abs(Z) < Inf)
  for(k in 1:K){
    # calculate/update the lower bound
    a[k] <- gsdmvn::gs_b(k = k, par = a #, par = lpar,
                         # info = info1, r = r, tol = tol, test_bound = test_lower, hgm1 = hgm1_1
                         # theta = rep(0, K), efficacy = FALSE)
    )
    # calculate/update the upper bound
    b[k] <- gsdmvn::gs_b(k = k, par = b #, par = upar,
                         #info = info0, r = r, tol = tol, test_bound = test_upper, hgm1 = hgm1_0,  
    )
    if(k == 1){
      # upperProb[1] <- if(b[1] < Inf) {pnorm(b[1], mean = sqrt(info[1]) * theta[1], lower.tail = FALSE)}else{0}
      # lowerProb[1] <- if(a[1] > -Inf){pnorm(a[1], mean = sqrt(info[1]) * theta[1])}else{0}
      upperProb[1] <- if(b[1] < Inf) {pnorm(b[1], mean = 0, lower.tail = FALSE)}else{0}
      lowerProb[1] <- if(a[1] > -Inf){pnorm(a[1], mean = 0)}else{0}
      # hgm1_1 <- h1(r = r, theta = theta1[1], I = info1[1], a = a[1], b = b[1])
      # hgm1   <- h1(r = r, theta = theta[1],  I = info[1],  a = a[1], b = b[1])
      hgm1   <- h1(r = r, theta = 0, I = info[1],  a = a[1], b = b[1])
    }else{
      # calculate the probability to cross upper bound
      upperProb[k] <- if(b[k] < Inf){
        # hupdate(r = r, theta = theta[k], I = info[k], a = b[k], b = Inf,
        #         thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
        hupdate(r = r, theta = 0, I = info[k], a = b[k], b = Inf,
                thetam1 = 0, Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
      }else{0}
      # calculate the probability to cross lower bound
      lowerProb[k] <- if(a[k] > -Inf){
        # hupdate(r = r, theta = theta[k], I = info[k], a = -Inf, b = a[k],
        #         thetam1 = theta[k - 1], Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
        hupdate(r = r, theta = 0, I = info[k], a = -Inf, b = a[k],
                thetam1 = 0, Im1 = info[k - 1], gm1 = hgm1) %>% summarise(sum(h)) %>% as.numeric()
      }else{0}
      
      if(k < K){
        # hgm1_0 <- hupdate(r = r, theta = 0, thetam1 = 0, I = info0[k], Im1 = info0[k-1],
        #                   a = if(binding){a[k]}else{-Inf}, b = b[k], gm1 = hgm1_0)
        # hgm1_1 <- hupdate(r = r, theta = theta1[k], thetam1 = theta1[k-1],I = info1[k], Im1 = info1[k-1],
        #                   a = a[k], b = b[k], gm1 = hgm1_1)
        # hgm1   <- hupdate(r = r, theta = theta[k], thetam1 = theta[k-1], I = info[k],  Im1 = info[k-1],
        #                   a = a[k], b = b[k], gm1 = hgm1)
        hgm1   <- hupdate(r = r, theta = 0, thetam1 = 0, I = info[k], Im1 = info[k-1],  
                          a = a[k], b = b[k], gm1 = hgm1)
      }
    }
  }
  
  table_H0 <- tibble::tibble(
    Analysis = rep(1:K, 2),
    Bound = c(rep("Upper", K), rep("Lower", K)),
    Z = c(b, a),
    Probability = c(cumsum(upperProb),
                    cumsum(lowerProb)),
    theta = rep(rep(0, K), 2),
    theta1 = rep(rep(0, K), 2),
    IF = rep(info1 / max(info1), 2),
    info = rep(info0, 2),  #rep(info, 2),
    info0 = rep(info0, 2),
    info1 = rep(info1, 2),
    hypothesis = rep("H0", 2*K))%>% filter(abs(Z) < Inf)
  
  out <- union_all(
    table_H1,
    table_H0
    ) %>% arrange(desc(hypothesis), desc(Bound), Analysis)
  return(out)
}