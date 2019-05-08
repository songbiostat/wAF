#' Score Statistics Permutation
#'
#' @description It generates correlated score statistics by permutation.
#'
#' @param Y Phenotype data. It can be a continuous trait
#' or a disease indicator. A vector or length n (number of subjects).
#' @param X Genotype data. A matrix with dimension n (number of subjects)
#' by K (number of variants).
#' @param binary Indicator of whether Y is binary.
#' @param out Use "2sided" for outputting P-values of two-sided score tests,
#' "all" for outputting P-values of both one-sided and two score tests.
#' @param nperm Number of permutations.
#'
#' @return A list object. Us: score statistics.
#' nonzero_var: indicators of variants with nonzero variance.
#' If out="2sided", pvs: P-values of two-tailed score tests.
#' If out="all", pvs_2sided: P-values of two_tailed score tests;
#' pvs_left: left-side P-values of one-tailed score tests.
#' pvs-right: right-side P-values of one-tailed score tests.
#' @export
#'
perm_score <- function(Y, X, binary = FALSE, out = c("2sided", "all"), nperm = 1e4) {

  out <- match.arg(out)

  Ybar <- mean(Y)
  resid <- Y - Ybar
  X.c <- scale(X, scale = FALSE)

  if (binary) {
    v <- diag(Ybar * (1 - Ybar) * t(X.c) %*% X.c)
  } else {
    v <- diag(var(Y) * t(X.c) %*% X.c)
  }

  resid.perm <- cbind(resid, replicate(nperm, sample(resid)))
  index <-v != 0
  Up <- ((t(X) %*% resid.perm)/sqrt(v))[index, ]
  p.2side <- 2 * (1 - pnorm(abs(Up)))

  if (out == "2sided") {
    result <- list(pvs = p.2side, Us = Up, nonzero_var = index)
  }else{
    p.right <- 1-pnorm(Up)
    p.left <- pnorm(Up)
    result <- list(pvs_2sided = p.2side, pvs_right = p.right, pvs_left = p.left,
                   Us = Up, nonzero_var = index)
  }
  return(result)
}
