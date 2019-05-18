#' Weighted Adaptive Fisher Test for Trait-SNV Set Association
#'
#' @description This function performs weighted Adaptive Fisher (wAF) test for
#' detecting association between a single trait and a set of single nucleotide
#' variatnts (SNVs).
#'
#' @param Y Phenotype data. It can be a continuous trait
#' or a binary trait. A vector or length n (number of subjects).
#'
#' @param X Genotype data. A matrix with dimensions n (number of subjects)
#' by K (number of variants).
#'
#' @param binary Indicator of whether Y is binary.
#'
#' @param method Method option. Use "wAF" for wAF test,
#' "wAFd" for directed wAF test.
#'
#' @param weight Weight option. Use "sd" for standard deviation weights,
#' "flat" for flat weights.
#'
#' @param weight_vector User-specified weights. A vector of length K (number of
#' variants).
#'
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#'
#'  @param nperm Number of permutations. Default is 10,000.
#'
#' @return An object of "wAF" class.
#' \describe{
#'  \item{pv}{P-value of wAF(wAFd) test.}
#'  \item{stat}{Test statistic of wAF(wAFd) test.}
#'  \item{SNV_combined}{Variants which are combined into the test
#'   statistic. The index of included variants are returned in the ascending
#'   order of their weighted P-values.}
#'   \item{method}{Method used.}
#'   \item{weight_method}{Method of weighing variants, "sd" of "flat".}
#'   \item{weight_vector}{Vector of Weights used (if "sd" or self-defined
#'   weights are used).}
#'
#' }
#'
#' @export
#'
#' @examples
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' test1 <- wAF(Y, X, method = "wAF", nperm = 100)
#' test2 <- wAF(Y, X, method = "wAFd", weight = "flat", nperm = 100)
#'
#' summary(test1)
#' test2
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
    weight <- "self defined"
  }

  if (method == "wAF") {
    test <- wAF_combine(score$pvs, weight = weight_vector, n0 = n0)
  } else {
    test <- wAFd_combine(score$pvs_2sided, score$pvs_right, score$pvs_left,
                         weight = weight_vector, n0 = n0)
  }

  combine.index <- ((1:ncol(X))[score$nonzero_var])[test$SNV_combined]

  result <- list(pv = test$pv, stat = test$stat, SNV_combined = combine.index,
                 method = method, weight_method = weight)

  if (weight %in% c("sd", "self defined")) {
    result[["weight_vector"]] <- weight_vector
  }

  class(result)<-"wAF"

  return(result)
}
