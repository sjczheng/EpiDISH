#' Dummy beta value matrix
#'
#' A dataset containing dummy whole blood beta value matrix
#' The dataset is part of the 450k GEO dataset GSE80559
#' The tissue type is whole blood
#' To reduce the data size, only 1000 probes are included
#' 330 probes are overlapped with centDHSblood.m
#' You can get the whole beta value matrix by
#' \code{exprs(GEOquery::getGEO('GSE80559')[[1]])}
#'
#' \itemize{
#'   \item beta value matrix of 1000 probes and 2 samples
#' }
#'
#' @docType data
#' @keywords datasets
#' @name DummyBeta.m
#' @usage data(DummyBeta.m)
#' @format A matrix with 1000 rows and 2 columns
NULL
