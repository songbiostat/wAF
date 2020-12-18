#' Score Statistics Permutation
#'
#' @description This function generates correlated score statistics
#' by permutation.
#'
#' @param Y Phenotype data. It can be a continuous trait
#' or a binary trait. A vector or length n (number of subjects).
#' @param X Genotype data. A matrix with dimensions n (number of subjects)
#' by K (number of variants). If elements within a column of X are identical,
#' this column needs to be removed ahead.
#' @param binary Indicator of whether Y is binary.
#' @param cov Covariates. A matrix with dimensions n (number of subjects)
#' by J (number of covariates).
#' @param nperm Number of permutations. Default is 1,000.
#' @param seed Specify seed for permutations.
#'
#' @return A list object.
#' \describe{
#'   \item{Us}{Score statistics. The first column contains the original
#'   data.}
#'   \item{pvs}{P-values of two-tailed score tests;
#'   can be further used by wAF_combine function.}
#'   \item{pvs_left}{Left-side P-values of one-tailed score tests;
#'   can be further used by wAF_combine, wAFd_combine functions.}
#'}
#'
#' @export
#'
#' @seealso {\code{\link{set.seed}}}
#'
#' @examples
#' # binary trait
#' Y <- RV_dense$trait
#' X <- RV_dense$SNV[, -RV_dense$zero_var]
#' result <- perm_score(Y, X, binary = TRUE, nperm = 100)
#' names(result)
#'
#' # continuous trait
#' Y <- SNV_sparse$trait
#' X <- SNV_sparse$SNV[, -SNV_sparse$zero_var]
#' result <- perm_score(Y, X, nperm = 100)
#' names(result)

perm_score <- function(Y, X, binary = FALSE, cov = NULL,
                       nperm = 1e3, seed = NULL) {
  K <- ncol(X)
  n <- nrow(X)

  if (length(Y) != n)
    stop("length of Y should be the same as row numbers of X")

  s <- apply(X, 2, sd)
  if (any (s == 0)) {
    ind <- which(s == 0)
    print("column indexes of X with no variation")
    print(ind)
    X <- X[, -ind]
  }

  if(is.null(cov)) cov <- rep(1, length(Y))

  ## fit the null model
  model <- ifelse(binary, "binomial", "gaussian")
  data1 <- data.frame(trait = Y, cov)
  null.model <- glm(trait ~ ., family = model, data = data1)

  cov <- data.frame(cov)
  X.null <- lm(X ~ ., data = cov)
  Xres <- resid(X.null)
  if(K == 1) {Xres <- matrix(Xres, ncol = K)}

  ## calculate score test statistics U
  U <- glm.scoretest(null.model, Xres)

  ## permute X residuals and calculate permutation Us
  if (!is.null(seed)) {
    set.seed(seed)
  }
  order.perm <- replicate(nperm, sample(1:n))
  U.perm <- apply(order.perm, 2,
                  function(x) glm.scoretest(null.model, Xres[x, ]))

  ## Calculate p-values
  Up <- cbind(U, U.perm)
  p <- 2 * (1 - pnorm(abs(Up)))
  p1 <- pnorm(Up)

  result <- list(Us = Up, pvs = p, pvs_left = p1)
  return(result)
  ## Us are the score statistics
  ## pvs are the p-values;
  ## pvs_left are left-sided P-values; can be further used by wAFd
}
