% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine.R
\name{wAFd_combine}
\alias{wAFd_combine}
\title{Directed wAF Combination}
\usage{
wAFd_combine(pvs_2sided, pvs_left, log = FALSE, weight = NULL, n0 = 1)
}
\arguments{
\item{pvs_2sided}{P-values of two-tailed tests. A matrix with
dimenstions K by N. If an object of P-values from perm_score is used,
K is the number of SNVs, and N is the number of permutations plus 1.}

\item{pvs_left}{Left-side P-values of one-tailed tests. Dimensions
should be the same as pvs_2sided.}

\item{log}{Indicator of whether P-values are on the log scale.}

\item{weight}{Weights given to the P-values. A vector with dimension
K. Flat weight is used if it is not specified.}

\item{n0}{Tuning parameter. Discard the first n0-1 P-values of each
column.}
}
\value{
A list object.
\describe{
  \item{pv}{P-value of wAF dtest.}
  \item{stat}{Observed test statitic of wAFd test.}
  \item{loci_combined}{Variants which are combined into the test
  statistic. The index of included variants are returned in the ascending
  order of their weighted P-values.}
  \item{stat_all}{wAFd statistics for all permuted samples.}
  \item{pv_all}{P-values of wAFd statistics for all permuted samples.}
}
}
\description{
This function combines P-values of one-tailed and
two-tailed tests using directed weighted adaptive Fisher (wAFd)
method.
}
\examples{
# Combine P-values of normally distributed test statistics
U <- matrix(rnorm(10 * 100), ncol=100)
p <- 2 * (1 - pnorm(abs(U)))
p.l <- pnorm(abs(U))
wt <- (1:10)/55
test <- wAFd_combine(p, p.l, weight = wt)

Y <- SNV_sparse$trait
X <- SNV_sparse$SNV[, - SNV_sparse$zero_var]
result <- perm_score(Y, X, nperm = 100)
test <- wAFd_combine(result$pvs, result$pvs_left)

}
