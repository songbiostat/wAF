#' Simulated Rare Variants Data in Dense Scenario
#'
#' A simulated dataset containing 1,000 subjects and 300 rare variants.
#' 20\% of the variants are simulated to be causal/associated. Effect
#' strengths are randomly sampled from U(-0.5, 0.5) distribution.
#'
#' @docType data
#'
#' @usage data(RV_dense)
#'
#' @format A list object.
#' \describe{
#'   \item{SNV}{A 1,000 by 300 matrix containing the genotypes. Each
#'   component of the matrix denotes the number of minor alleles.}
#'   \item{trait}{A vector of length 1,000 containing disease labels
#'   (500 cases and 500 controls).}
#'   \item{zero_var}{Indexes for columns with no variation. These
#'   columns should be removed if SNV is further used by perm_score, wAF and
#'   wAFd functions.}
#' }
"RV_dense"

#' Simulated Rare Variants Data in Sparse Scenario
#'
#' A simulated dataset containing 1,000 subjects and 300 rare variants.
#' 1\% of the variants are simulated to be causal/associated. Effect
#' strengths are randomly sampled from U(-2, 2) distribution.
#'
#' @docType data
#'
#' @usage data(RV_sparse)
#'
#' @format A list object.
#' \describe{
#'   \item{SNV}{A 1,000 by 300 matrix containing the genotypes. Each
#'   component of the matrix denotes the number of minor alleles.}
#'   \item{trait}{A vector of length 1,000 containing disease labels
#'   (500 cases and 500 controls).}
#'   \item{zero_var}{Indexes for columns with no variation. These
#'   columns should be removed if SNV is further used by perm_score, wAF and
#'   wAFd functions.}
#' }
"RV_sparse"

#' Simulated Single Nucleotide Variants (SNVs) Data in Dense Scenario
#'
#' A simulated dataset containing 1,000 subjects and 100 SNVs. 20\% of
#' the variants are simulated to be causal/associated. Effect strengths
#' are randomly sampled from U(-0.2, 0.2) distribution.
#'
#' @docType data
#'
#' @usage data(SNV_dense)
#'
#' @format A list object.
#' \describe{
#'   \item{SNV}{A 1,000 by 100 matrix containing the genotypes. Each
#'   component of the matrix denotes the number of minor alleles.}
#'   \item{trait}{A vector of length 1,000 containing a continuous trait.}
#'   \item{zero_var}{Indexes for columns with no variation. These
#'   columns should be removed if SNV is further used by perm_score, wAF and
#'   wAFd functions.}
#' }
"SNV_dense"

#' Simulated Single Nucleotide Variants (SNVs) Data in Sparse Scenario
#'
#' A simulated dataset containing 1,000 subjects and 100 SNVs. 1\% of
#' the variants are simulated to be causal/associated. Effect strengths
#' are randomly sampled from U(-0.5, 0.5) distribution.
#'
#' @docType data
#'
#' @usage data(SNV_sparse)
#'
#' @format A list object.
#' \describe{
#'   \item{SNV}{A 1,000 by 300 matrix containing the genotypes. Each
#'   component of the matrix denotes the number of minor alleles.}
#'   \item{trait}{A vector of length 1,000 containing a continuous trait.}
#'   \item{zero_var}{Indexes for columns with no variation. These
#'   columns should be removed if SNV is further used by perm_score, wAF and
#'   wAFd functions.}
#' }
"SNV_sparse"
