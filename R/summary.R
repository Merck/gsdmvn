summary <- function(x, ...) {
  UseMethod("summary", x)
}

#' Title
#' @title summary for \code{fixed_design()} object
#' @param x a fixed design object returned by \code{fixed_design()}
#' @param ... additional arguments
#' @export summary.fixed_design
#' @exportS3Method 
#' @examples

summary.fixed_design <- function(x, ...){
  x_design <- switch(x$design,
                     "AHR" = {"AHR"},
                     "LF" = {"LF"},
                     "RD" = {"RD"},
                     "FH" = {paste0("FH with (rho = ", x$design_par$rho, 
                                    ", gamma = ", gamma = x$design_par$gamma, ")")},
                     "MB" = {paste0("MB with (rho = ", x$design_par$rho, 
                                    ", gamma = ",  x$design_par$gamma, 
                                    ", tau =  ", x$design_par$tau, ")")},  
                     "MaxCombo" = {paste0("MaxCombo with (rho, gamma, tau) = (",
                                          paste(apply(do.call(rbind, x$design_par), 2 , paste , collapse = "," ), collapse = ") and ("),
                                          ")")}
                     )
  
  ans <- x$analysis %>% mutate(Design = x_design)
  class(ans) <- c("fixed_design", x$design, class(ans))
  return(ans)
}
