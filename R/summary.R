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
#' X <- SNV_sparse$SNV[, -SNV_sparse$zero_var]
#' test1 <- wAF(Y, X, nperm = 100)
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
  if (x$weight == "flat") {
    cat("Flat Weights\n")
  } else {
    cat(paste(x$weight, "weights: \n"))
    print(x$weight_values)
  }
  cat("\n")
  cat("SNVs combined into test statistic:\n")
  print(x$loci_combined)
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
#' X <- SNV_sparse$SNV[, -SNV_sparse$zero_var]
#' test2 <- wAF(Y, X, w = "flat", nperm = 100)
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
  if (x$weight == "flat") {
    cat("Flat Weights\n")
  } else {
    cat(paste(x$weight, "weights: \n"))
    print(x$weight_values)
  }
  cat("\n")
  cat("SNVs combined into test statistic:\n")
  print(x$loci_combined)
}
