\name{pantissue8CTmsa.m}
\alias{pantissue8CTmsa.m}
\docType{data}
\title{Pan-tissue DNA methylation reference for eight cell types on MSA}
\description{
A reference matrix for estimating cell-type fractions with \code{\link{epidish}}.
This reference uses MSA probe identifiers with suffixes.
Match the row names of the sample matrix to the identifiers used in this reference.
}
\usage{data(pantissue8CTmsa.m)}
\format{
An 800 by 8 numeric matrix of DNA methylation beta values. The columns
are Epithelial, Endothelial, Fibroblast, Adipocyte, B, non-B lymphoid,
MonoMac and Granulocyte.
}
\seealso{\code{\link{pantissue8CT.m}}, \code{\link{pantissue8CTepic2.m}}}
\keyword{datasets}
