#' Score Statistics Permutation
#'
#' @description This function generates correlated score statistics
#' by permutation.
#'
#' @param Y Phenotype data. It can be a continuous trait
#' or a binary trait. A vector or length n (number of subjects).
#'
#' @param X Genotype data. A matrix with dimensions n (number of subjects)
#' by K (number of variants).
#'
#' @param binary Indicator of whether Y is binary.
#'
#' @param out Output options. "2sided" only outputs P-values of two-tailed
#' score tests; "all" outputs P-values of two-tailed test and two
#' one-tailed score tests.
#'
#' @param nperm Number of permutations. Default is 1,000.
#'
#' @return A list object.
#' \describe{
#'   \item{Us}{Score statistics. The first column contains the original
#'   data.}
#'   \item{nonzero_var}{Indicators of variants with nonzero variance.}
#'   \item{pvs}{P-values of two-tailed score tests (if out = "2sided" is
#' used).}
#'   \item{pvs_2sided}{P-values of two_tailed score tests (if out = "all"
#'   is used).}
#'   \item{pvs_left}{Left-side P-values of one-tailed score tests
#' (if out = "all" is used).}
#'   \item{pvs_sided}{Right-side P-values of one-tailed score tests
#' (if out = "all" is used).}
#'}
#'
#' @export
#'
#' @examples
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV
#' result <- perm_score(Y, X, out = "all", nperm = 100)
#' names(result)
perm_score <- function(Y, X, binary = FALSE,
                       out = c("2sided", "all"), nperm = 1e3) {

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
    result <- list(pvs_2sided = p.2side, pvs_right = p.right,
                   pvs_left = p.left, Us = Up, nonzero_var = index)
  }
  return(result)
}
