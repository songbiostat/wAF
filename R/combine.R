#' wAF Combination
#'
#' @description This function combines P-values using
#' weighted adaptive Fisher (wAF) method.
#'
#' @param p P-values to be combined. A matrix with dimenstions K by N.
#' If an object of P-values from perm_score is used, K is the number
#' of SNVs, and N is the number of permutations plus 1.
#'
#' @param log Indicator of whether P-values are on the log scale.
#'
#' @param weight Weights given to the P-values. A vector with dimension
#' K. Flat weight is used if it is not specified.
#'
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each
#' column.
#'
#' @return A list object.
#' \describe{
#'   \item{pv}{P-value of wAF test.}
#'   \item{stat}{Observed test statitic of wAF test.}
#'   \item{stat_dist}{wAF test statistics for P-values in every column.}
#'   \item{SNV_combined}{Variants which are combined into the test
#'   statistic. The index of included variants are returned in the ascending
#'   order of their weighted P-values.}
#' }
#'
#' @export
#'
#' @examples
#' # Combine P-values of normally distributed test statistics
#' U <- matrix(rnorm(10 * 100), ncol=100)
#' p <- 2 * (1 - pnorm(abs(U)))
#' wt <- (1:10)/55
#' test <- wAF_combine(p, weight = wt)
#'
#'
#' # Combine P-values from perm_score
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' result <- perm_score(Y, X, out = "all", nperm = 100)
#' test <- wAF_combine(result$pvs_2sided)
#'
wAF_combine <- function(p, log = FALSE, weight = NULL, n0 = 1) {

  N <- nrow(p)
  T <- ncol(p)

  if (is.null(weight)) {weight <- rep(1, N)}

  if (log == FALSE) {
    r <- weight * log(p)
  } else {
    r <- weight * p
  }

  s <- -apply(apply(r, 2, sort), 2, cumsum)
  p.perm <- t(apply(-s, 1, rank, ties.method="max")/T)
  stat <- apply(p.perm[n0:N,], 2, min)
  p.wAF <- mean(stat <= stat[1])

  # SNVs combined in wAF statistic
  obs.order <- order(r[,1])
  combine.index <- obs.order[1:which.min(p.perm[n0:N, 1])]

  result <- list(pv = p.wAF, stat = stat[1], stat_dist = stat,
                 SNV_combined = combine.index)

  return(result)
}


#' Directed wAF Combination
#'
#' @description This function combines P-values of one-tailed and
#' two-tailed tests using directed weighted adaptive Fisher (wAFd)
#' method.
#'
#' @param pvs_2sided P-values of two-tailed tests. A matrix with
#' dimenstions K by N. If an object of P-values from perm_score is used,
#' K is the number of SNVs, and N is the number of permutations plus 1.
#'
#' @param pvs_right Right-side P-values of one-tailed tests. Dimensions
#' should be the same as pvs_2sided.
#'
#' @param pvs_left Left-side P-values of one-tailed tests. Dimensions
#' should be the same as pvs_2sided.
#'
#' @param log Indicator of whether P-values are on the log scale.
#'
#' @param weight Weights given to the P-values. A vector with dimension
#' K. Flat weight is used if it is not specified.
#'
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each
#' column.
#'
#' @return A list object.
#' \describe{
#'   \item{pv}{P-value of wAF dtest.}
#'   \item{stat}{Observed test statitic of wAFd test.}
#'   \item{stat_dist}{wAFd test statistics for P-values in every column.}
#'   \item{SNV_combined}{Variants which are combined into the test
#'   statistic. The index of included variants are returned in the ascending
#'   order of their weighted P-values.}
#' }
#'
#' @export
#'
#' @examples
#' # Combine P-values of normally distributed test statistics
#' U <- matrix(rnorm(10 * 100), ncol=100)
#' p <- 2 * (1 - pnorm(abs(U)))
#' p.r <- 1 - pnorm(abs(U))
#' p.l <- pnorm(abs(U))
#' wt <- (1:10)/55
#' test <- wAFd_combine(p, p.r, p.l, weight = wt)
#'
#'
#' # Combine P-values from perm_score
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' result <- perm_score(Y, X, out = "all", nperm = 100)
#' test <- wAFd_combine(result$pvs_2sided, result$pvs_right,
#' result$pvs_left)
#'
wAFd_combine <- function(pvs_2sided, pvs_right, pvs_left,
                         log = FALSE, weight = NULL, n0 = 1){
  wAF.2 <- wAF_combine(pvs_2sided, log = log,
                       weight = weight, n0 = n0)
  wAF.r <- wAF_combine(pvs_right, log = log,
                       weight = weight, n0 = n0)
  wAF.l <- wAF_combine(pvs_left, log = log,
                       weight = weight, n0 = n0)

  stat.2side <- wAF.2$stat_dist
  stat.right <- wAF.r$stat_dist
  stat.left <- wAF.l$stat_dist

  stat.all <- rbind(stat.2side, stat.right, stat.left)
  min.all <- apply(stat.all, 2, min)
  p.wAFd <- mean(min.all <= min.all[1])

  # SNVs combined in wAF statistic
  combine.list <- list(wAF.2$SNV_combined, wAF.r$SNV_combined,
                       wAF.l$SNV_combined)
  combine.index <- combine.list[[which.min(stat.all[,1])]]

  result <- list(pv = p.wAFd, stat = min.all[1], stat_dist = min.all,
                 SNV_combined = combine.index)

  return(result)
}
