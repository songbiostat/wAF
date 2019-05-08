#' Weighted Adaptive Fisher Test for A Single Trait-SNV Set Association
#'
#' @description It performs weighted adaptive Fisher for a single trait and
#' a set of SNV data.
#'
#' @param Y Phenotype data. It can be a continuous trait
#' or a disease indicator. A vector or length n (number of subjects).
#' @param X X Genotype data. A matrix with dimension n (number of subjects)
#' by K (number of variants).
#' @param binary Indicator of whether Y is binary.
#' @param method Use "wAF" for weighted adaptive Fisher test, "wAFd" for directed
#' adaptive Fisher test.
#' @param weight Use "sd" for standard deviation weights, "flat" for flat weights.
#' @param weight_vector User-specified weights. A vector of length K (number of
#' variants).
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#' @param nperm Number of permutations.
#'
#' @return A list object. pv: P-value of wAF(wAFd) test;
#' stat: observed test statitic of wAF(wAFd) test;
#' stat_dist: wAF(wAFd) test statistics for P-values in every column.
#' @export
#'
wAF <- function(Y, X, binary = FALSE, method = c("wAF", "wAFd"),
                weight = c("sd", "flat"), weight_vector = NULL,
                n0 = 1, nperm = 1e4 ){
  method <- match.arg(method)
  weight <- match.arg(weight)

  out <- ifelse(method == "wAF", "2sided", "all")

  score <- perm_score(Y, X, binary = binary, out = out, nperm = nperm)

  if (is.null(weight_vector)) {
    if (weight == "sd") {
      weight_vector <- (colMeans(X)/2)[score$nonzero_var]
    }
  } else {
    weight_vector <- weight_vector[score$nonzero_var]
  }

  if (method == "wAF") {
    test <- wAF_combine(score$pvs, weight = weight_vector, n0 = n0)
  } else {
    test <- wAFd_combine(score$pvs_2sided, score$pvs_right, score$pvs_left,
                         weight = weight_vector, n0 = n0)
  }

  result <- list(pv = test$pv, stat = test$stat, stat_dist = test$stat_dist)
  return(result)
}
