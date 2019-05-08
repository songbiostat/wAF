#' wAF Combination
#'
#' @description It combines P-values using weighted adaptive Fisher (wAF) method.
#'
#' @param p P-values to combine. A matrix with dimenstion K (number of variants)
#' by N (number of permutations plus one).
#' @param log Indicator of whether P-values are on the log scale.
#' @param weight Weights given to the P-values. A vector with dimention
#' K (number of variants).Flat weight is used if it is not specified.
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#'
#' @return A list object. pv: P-value of wAF test; stat: observed test statitic of wAF
#' test; stat_dist: wAF test statistics for P-values in every column.
#' @export
#'
wAF_combine <- function(p, log = FALSE, weight = NULL, n0 = 1) {

  N <- nrow(p)
  T <- ncol(p)

  if (is.null(weight)) {weight <- rep(1, N)}

  if (log == FALSE) {
    r <- weight*log(p)
  } else {
    r <- weight*p
  }

  s <- -apply(apply(r,2,sort),2,cumsum)
  p.perm <- t(apply(-s, 1, rank, ties.method="max")/T)
  stat <- apply(p.perm[n0:N,], 2, min)
  p.wAF <- mean(stat <= stat[1])
  result <- list(pv = p.wAF, stat = stat[1], stat_dist = stat)

  return(result)
}


#' wAFd Combination
#'
#' @description It combines P-values of one-tailed and two-tailed tests
#' using directed weighted adaptive Fisher (wAFd) method.
#' @param pvs_2sided P-values of two-tailed tests. A matrix with dimenstion K (number of variants)
#' by N (number of permutations plus one).
#' @param pvs_right Right-side P-values of one-tailed tests. A matrix with dimenstion K (number of variants)
#' by N (number of permutations plus one).
#' @param pvs_left Left-side P-values of one-tailed tests. A matrix with dimenstion K (number of variants)
#' by N (number of permutations plus one).
#' @param log Indicator of whether P-values are on the log scale.
#' @param weight Weights given to the P-values. A vector with dimention
#' K (number of variants).Flat weight is used if it is not specified.
#' @param n0 Tuning parameter. Discard the first n0-1 P-values of each column.
#'
#' @return A list object. pv: P-value of wAFd test; stat: observed test statitic of wAFd
#' test; stat_dist: wAFd test statistics for P-values in every column.
#' @export
#'
wAFd_combine <- function(pvs_2sided, pvs_right, pvs_left,
                         log = FALSE, weight = NULL, n0 = 1){
  stat.2side <- wAF_combine(pvs_2sided, log = log,
                            weight = weight, n0 = n0)$stat_dist
  stat.right <- wAF_combine(pvs_right, log = log,
                            weight = weight, n0 = n0)$stat_dist
  stat.left <- wAF_combine(pvs_left, log = log,
                           weight = weight, n0 = n0)$stat_dist

  stat.all <- rbind(stat.2side,stat.right,stat.left)
  min.all <- apply(stat.all,2,min)
  p.wAFd <- mean(min.all<=min.all[1])
  result <- list(pv = p.wAFd, stat = min.all[1], stat_dist = min.all
                 , stat.2side,stat.right,stat.left)
  # comment the previous line for the final version

  return(result)
}
