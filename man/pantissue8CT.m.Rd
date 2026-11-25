\name{pantissue8CT.m}
\alias{pantissue8CT.m}
\docType{data}
\title{Pan-tissue DNA methylation reference for eight cell types on 450K and EPIC v1}
\description{
A reference matrix for estimating cell-type fractions with \code{\link{epidish}}.
This reference is intended for Illumina 450K and EPIC v1 arrays.
Row names use unsuffixed CpG probe identifiers.
Match the row names of the sample matrix to the identifiers used in this reference.
}
\usage{data(pantissue8CT.m)}
\format{
An 800 by 8 numeric matrix of DNA methylation beta values. The columns
are Epithelial, Endothelial, Fibroblast, Adipocyte, B, non-B lymphoid,
MonoMac and Granulocyte.
}
\seealso{\code{\link{pantissue8CTepic2.m}}, \code{\link{pantissue8CTmsa.m}}}
\keyword{datasets}
