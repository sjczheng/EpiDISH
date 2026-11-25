\name{calculate_mFDI}
\alias{calculate_mFDI}
\title{Calculate the extrinsic mFDI score}
\usage{calculate_mFDI(estF, epsilon = 0.01)}
\arguments{
  \item{estF}{An \code{\link{epidish}} result or a matrix or data frame
  of numeric fractions with samples in rows and columns named
  \code{Fibroblast} and \code{Epithelial}.}
  \item{epsilon}{A numeric pseudocount added to both fractions.}
}
\description{
Calculates the base-2 logarithm of the ratio of the fibroblast fraction
plus the pseudocount to the epithelial fraction plus the pseudocount.
}
\value{
A numeric vector named by sample row names when available.
Entries with non-finite fractions or nonpositive adjusted fractions are
returned as \code{NA}.
}
\examples{
fractions <- cbind(Fibroblast = c(0.2, 0.6), Epithelial = c(0.6, 0.2))
calculate_mFDI(fractions)
}
