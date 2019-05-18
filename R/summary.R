#' Summary Function for Objects of "wAF" Class
#'
#' @param x An object of "wAF" class.
#'
#' @return Method used; P-value; weights used; variants combined into the test
#' statistic.
#'
#' @export
#'
#' @examples
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' test1 <- wAF(Y, X, method = "wAF", nperm = 100)
#' summary(test1)
#'
summary.wAF <- function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Weights:\n")
  if (x$weight_method == "flat") {
    cat("Flat Weights\n")
  } else {
    cat(paste(x$weight_method, "weights: \n"))
    print(x$weight_vector)
  }
  cat("\n")
  cat("SNVs combined into test statistic:\n")
  print(x$SNV_combined)
}


#' Print Function for Objects of "wAF" Class
#'
#' @param x An object of "wAF" class.
#'
#' @return Method used; P-value; weights used; variants combined into the test
#' statistic.
#'
#' @export
#'
#' @examples
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' test2 <- wAF(Y, X, method = "wAFd", weight = "flat", nperm = 100)
#' test2
#' print(test2)
#'
print.wAF <- function(x, ...){
  cat("Method:\n")
  cat(paste(x$method, "\n"))
  cat("\n")
  cat("P-value:\n")
  print(x$pv)
  cat("\n")
  cat("Weights:\n")
  if (x$weight_method == "flat") {
    cat("Flat Weights\n")
  } else {
    cat(paste(x$weight_method, "weights: \n"))
    print(x$weight_vector)
  }
  cat("\n")
  cat("SNVs combined into test statistic:\n")
  print(x$SNV_combined)
}
