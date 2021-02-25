#' P Value Flag Function
#'
#' This function is used to flag p values.
#' @param x P value.
#' @export


star.pval <- function(x, digits = 6, ...) {
  if (x > 0.001){
    sig_level <- cut(x,
                     breaks = c(0.001, 0.01, 0.05, Inf),
                     labels = c("**", "*", ""))
    # Next, we print the values, keeping only the set number of significant digits,
    # appending the asterisks, and making sure we don't include quotes so they
    # still look like numbers.
    out <- paste0(signif(x, digits = digits), sig_level)
    return(out)
  } else {
    return(paste('<', 0.001, '***', sep = ''))
  }

}

