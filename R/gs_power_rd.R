gs_power_rd <- function(){
  x <- gs_info_ahr(enrollRates = enrollRates,
                   failRates = failRates,
                   ratio = ratio,
                   events = events,
                   analysisTimes = analysisTimes
  )
  return(gs_power_npe(theta = x$theta, info = x$info, info0 = x$info0, binding = binding,
                      upper=upper, lower=lower, upar = upar, lpar= lpar,
                      test_upper = test_upper, test_lower = test_lower,
                      r = r, tol = tol) %>%
           right_join(x %>% select(-c(info, info0, theta)), by = "Analysis") %>%
           select(c(Analysis, Bound, Time, Events, Z, Probability, AHR, theta, info, info0)) %>%
           arrange(desc(Bound), Analysis)
  )
}