#' Directed Weighted Adaptive Fisher Test for Trait-SNV Set Association
#'
#' @description This function performs weighted Adaptive Fisher (wAF) test for
#' detecting association between a single trait and a set of single nucleotide
#' variatnts (SNVs).
#'
#' @param Y Y Phenotype data. It can be a continuous trait
#' or a binary trait. A vector or length n (number of subjects).
#' @param X Genotype data. A matrix with dimensions n (number of subjects)
#' by K (number of variants).
#' @param binary Indicator of whether Y is binary.
#' @param cov Covariates. A matrix with dimensions n (number of subjects)
#' by J (number of covariates).
#' @param w Weight option. Use "sd" for standard deviation weights,
#' "flat" for flat weights.
#' @param weight User-specified weights. A vector of length K (number of
#' variants).
#' @param adapt_perm Whether "step-up" algorithm is used for P-value
#' calculation. If FALSE, function permutes nperm times and stops.
#' If TRUE, nperm will be increased 10 times each round if P-value
#' <= 5/nperm. Algorithm stops if P-value > 5/nperm or <= cutoff.
#' @param cutoff Cutoff for "step-up" algorithm.
#' @param nperm Number of permutations. Also the starting number of
#' permutations for "step-up" algorithm. Default is 1,000.
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#' @param seed Specify seed for permutations.
#'
#' @return An object of "wAF" class.
#' \describe{
#'  \item{pv}{P-value of wAF test.}
#'  \item{stat}{Test statistic of wAF test.}
#'  \item{loci_combined}{Variants which are combined into the test
#'   statistic. The index of included variants are returned in the ascending
#'   order of their weighted P-values.}
#'   \item{stat_all}{wAFd statistics for all permuted samples.}
#'   \item{pv_all}{P-values of wAFd statistics for all permuted samples.}
#'   \item{method}{Method used.}
#'   \item{weight}{Method of weighing variants, "sd" of "flat".}
#'   \item{weight_values}{Vector of weights used (if "sd" or user-specified
#'   weights are used).}
#' }
#'
#' @export
#'
#' @seealso \code{\link{set.seed}}
#'
#' @examples
#' ## Binary trait
#' Y <- RV_sparse$trait
#' X <- RV_sparse$SNV[, -RV_sparse$zero_var]
#'
#' # sd weights
#' test1 <- wAFd(Y, X, binary = TRUE, nperm = 100)
#' summary(test1)
#'
#' # flat weights
#' test2 <- wAFd(Y, X, w = "flat", nperm = 100, adapt_perm = TRUE)
#' test2
#'
wAFd <- function(Y, X, binary = FALSE, cov = NULL,
                 w = c("sd", "flat"), weight = NULL,
                 adapt_perm = FALSE, cutoff = 2.5e-6, nperm = 1e3,
                 n0 = 1, seed = NULL){

  w <- match.arg(w)
  if (is.null(weight)) {
    if (w == "sd") {
      maf <- colMeans(X)/2
      weight <- sqrt(maf * (1-maf))
    }
  } else {
    w <- "user-specified"
  }

  score <- perm_score(Y, X, binary = binary, cov = cov,
                      nperm = nperm, seed = seed)
  test <- wAFd_combine(score$pvs, score$pvs_left, weight = weight, n0 = n0)

  if (adapt_perm) {
    pvs2.all <- score$pvs
    pvsl.all <- score$pvs_left
    while(test$pv <= 5/nperm & test$pv > cutoff) {
      nperm.add <- nperm * 9
      nperm <- nperm + nperm.add
      print(paste("start", nperm, "permutations"))
      score.add <- perm_score(Y, X, binary = binary, cov = cov,
                              nperm = nperm.add)
      pvs2.all <- cbind(pvs2.all, score.add$pvs[,-1])
      pvsl.all <- cbind(pvsl.all, score.add$pvs_left[,-1])
      test <- wAFd_combine(pvs2.all, pvsl.all, weight = weight, n0 = n0)
    }
  }

  result <- list(pv = test$pv, stat = test$stat,
                 loci_combined = test$loci_combined,
                 pv_all = test$pv_all, stat_all = test$stat_all,
                 method = "wAFd", weight = w)

  if (w %in% c("sd", "user-specified")) {
    result[["weight_values"]] <- weight
  }

  class(result)<-"wAF"

  return(result)
}
