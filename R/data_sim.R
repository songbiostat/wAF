#' Data Simulation
#'
#' @param n Number of subjects.
#' @param K Number of SNVs.
#' @param prop Proportion of true causal/associated SNVs.
#' @param strg Effect size. SNV effects are selected from U(-strg, strg)
#' distribution.
#' @param binary Idicator of whether the trait is binary.
#' @param case Number of cases if binary trait is simulated. Default is n/2.
#' @param corr Correlation matrix of the underlying multivariate normal
#' distribution. Default is a first order autoregressive structure with
#' adjacent correlation 0.9.
#' @param maf Minor allele frequencies for the SNVs.
#' @param beta SNV effects including the constant. If beta is given, K,
#' prop and strg should be missing; otherwise, beta will be overwritten.
#'
#' @return A list object. SNV is a n by K matrix with counts of minor alleles,
#' trait is a n-vector with the simulated traits.
#' @export
#'
data_sim <- function(n, K, prop, strg, binary = FALSE, case = NULL,
                     corr = NULL, maf = NULL, beta = NULL){

  mis <- missing(K) + missing(prop) + missing(strg)

  if (mis == 0) {
    beta <- c(log(0.05/0.95), rep(0, K))
    K.sig <- floor(K * prop)
    sig.ind <- sample(2:(K+1), K.sig)
    beta[sig.ind] <- runif(K.sig, -strg, strg)
  } else {
    K <- length(beta)-1
  }

  if (is.null(corr)) {
    corr <- 0.9 ^ abs(outer(1:K, 1:K, "-"))
  }
  corr.sqrt <- t(chol(corr))

  if (is.null(maf)) {
    maf <- exp(runif(K, log(0.001), log(0.5)))
  }
  q <- qnorm(maf)


  Y<-vector() # trait
  X<-vector() # SNV
  if (binary == FALSE) {
    for (l in 1:n) {
      z1 <- t(corr.sqrt) %*% rnorm(K)
      z2 <- t(corr.sqrt) %*% rnorm(K)
      x <- (z1 < q) + (z2 < q)
      y <- sum(beta * c(1,x)) + rnorm(1)
      X <- cbind(X,x)
      Y <- c(Y,y)
    }
  } else {
    if (is.null(case)) case <- floor(n/2)
    start <- rep(0, 2) # count controls and cases
    end <- c(n-case, case) # numbers need to be reached
    while (sum(start < end) > 0) {
      z1 <- t(corr.sqrt) %*% rnorm(K)
      z2 <- t(corr.sqrt) %*% rnorm(K)
      x <- (z1 < q) + (z2 < q)
      y <- rbinom(1, 1, inv.logit(sum(beta * c(1,x))))
      if (start[1] >= end[1]) {
        while(y == 0) {
          z1 <- t(corr.sqrt) %*% rnorm(K)
          z2 <- t(corr.sqrt) %*% rnorm(K)
          x<-(z1 < q) + (z2 < q)
          y<-rbinom(1, 1, inv.logit(sum(beta * c(1,x))))
        }
      } else if(start[2] >= end [2]){
        while(y == 1) {
          z1 <- t(corr.sqrt) %*% rnorm(K)
          z2 <- t(corr.sqrt) %*% rnorm(K)
          x <- (z1 < q) + (z2 < q)
          y <- rbinom(1, 1, inv.logit(sum(beta * c(1,x))))
        }
      }
      X <- cbind(X,x)
      Y <- c(Y,y)
      start[y+1] <- start[y+1] + 1
    }
  }
  X<-t(X)

  result <- list(SNV = X, trait = Y)
  return (result)
}
