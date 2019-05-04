#' Whole blood reference of 188 tsDHS-DMCs and 7 blood cell subtypes
#'
#' This reference is a subset of \code{centDHSbloodDMC.m}, and contains 188 
#' DMCs which exhibit similar median DNAm values across epithelial cells, 
#' fibroblasts and ICs to ensure that the estimation of IC subtype fractions is 
#' not confounded by the epithelial and fibroblast cells in the sample. It should 
#' be used in the \code{hepidish} function to estimate fractions of immunce cell 
#'  subtypes.
#'
#' \itemize{
#'   \item B-cells
#'   \item CD4+ T-cells
#'   \item CD8+ T-cells 
#'   \item NK-cells
#'   \item Monocytes
#'   \item Neutrophils 
#'   \item Eosinophils
#' }
#'
#' @docType data
#' @keywords datasets
#' @name centBloodSub.m
#' @usage data(centBloodSub.m)
#' @format A matrix with 188 rows and 7 columns
#' @references 
#' Zheng SC, Webster AP, Dong D, Feber A, Graham DG, Sullivan R, Jevons S, Lovat LB, 
#' Beck S, Widschwendter M, Teschendorff AE
#' \emph{A novel cell-type deconvolution algorithm reveals substantial contamination by immune cells in saliva, buccal and cervix.}
#' Epigenomics (2018) 10: 925-940.
#' doi:\href{https://doi.org/10.2217/epi-2018-0037}{
#' 10.2217/epi-2018-0037}.
#' 
#' Teschendorff AE, Breeze CE, Zheng SC, Beck S. 
#' \emph{A comparison of reference-based algorithms for correcting cell-type 
#' heterogeneity in Epigenome-Wide Association Studies.}
#' BMC Bioinformatics (2017) 18: 105.
#' doi:\href{https://doi.org/10.1186/s12859-017-1511-5}{
#' 10.1186/s12859-017-1511-5}.
#' 
#' Reinius LE, Acevedo N, Joerink M, Pershagen G, Dahlen S-E, Greco D, 
#' Soderhall C, Scheynius A, Kere J.
#' \emph{Differential DNA methylation in purified human blood cells: 
#' implications for cell lineage and studies on disease susceptibility.}
#' PLoS ONE (2012) 7: e41361.
#' doi:\href{https://doi.org/10.1371/journal.pone.0041361}{
#' 10.1371/journal.pone.0041361}.
#' 
NULL
