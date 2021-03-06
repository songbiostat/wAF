% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wAFd.R
\name{wAFd}
\alias{wAFd}
\title{Directed Weighted Adaptive Fisher Test for Trait-SNV Set Association}
\usage{
wAFd(
  Y,
  X,
  binary = FALSE,
  cov = NULL,
  w = c("sd", "flat"),
  weight = NULL,
  adapt_perm = FALSE,
  cutoff = 2.5e-06,
  nperm = 1000,
  n0 = 1,
  seed = NULL
)
}
\arguments{
\item{Y}{Y Phenotype data. It can be a continuous trait
or a binary trait. A vector or length n (number of subjects).}

\item{X}{Genotype data. A matrix with dimensions n (number of subjects)
by K (number of variants).}

\item{binary}{Indicator of whether Y is binary.}

\item{cov}{Covariates. A matrix with dimensions n (number of subjects)
by J (number of covariates).}

\item{w}{Weight option. Use "sd" for standard deviation weights,
"flat" for flat weights.}

\item{weight}{User-specified weights. A vector of length K (number of
variants).}

\item{adapt_perm}{Whether "step-up" algorithm is used for P-value
calculation. If FALSE, function permutes nperm times and stops.
If TRUE, nperm will be increased 10 times each round if P-value
<= 5/nperm. Algorithm stops if P-value > 5/nperm or <= cutoff.}

\item{cutoff}{Cutoff for "step-up" algorithm.}

\item{nperm}{Number of permutations. Also the starting number of
permutations for "step-up" algorithm. Default is 1,000.}

\item{n0}{Tuning parameter. Discard the first n0-1 P-values of each column.}

\item{seed}{Specify seed for permutations.}
}
\value{
An object of "wAF" class.
\describe{
 \item{pv}{P-value of wAF test.}
 \item{stat}{Test statistic of wAF test.}
 \item{loci_combined}{Variants which are combined into the test
  statistic. The index of included variants are returned in the ascending
  order of their weighted P-values.}
  \item{stat_all}{wAFd statistics for all permuted samples.}
  \item{pv_all}{P-values of wAFd statistics for all permuted samples.}
  \item{method}{Method used.}
  \item{weight}{Method of weighing variants, "sd" of "flat".}
  \item{weight_values}{Vector of weights used (if "sd" or user-specified
  weights are used).}
}
}
\description{
This function performs weighted Adaptive Fisher (wAF) test for
detecting association between a single trait and a set of single nucleotide
variatnts (SNVs).
}
\examples{
## Binary trait
Y <- RV_sparse$trait
X <- RV_sparse$SNV[, -RV_sparse$zero_var]

# sd weights
test1 <- wAFd(Y, X, binary = TRUE, nperm = 100)
summary(test1)

# flat weights
test2 <- wAFd(Y, X, w = "flat", nperm = 100, adapt_perm = TRUE)
test2

}
\seealso{
\code{\link{set.seed}}
}
