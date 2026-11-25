\name{calculate_EMTS}
\alias{calculate_EMTS}
\title{Calculate the cell-intrinsic EpiEMT score}
\usage{calculate_EMTS(beta)}
\arguments{
  \item{beta}{A numeric matrix of DNA methylation beta values with
  unsuffixed CpG identifiers as row names and samples in columns.}
}
\description{
Standardizes each available CpG from a fixed 33-CpG signature across the
supplied cohort, applies its fixed age-effect direction, and averages the
resulting signed z-scores for each sample.
}
\details{
Non-finite values are omitted. CpGs with fewer than two finite observations
or zero variance are excluded. Scores depend on the cohort used for
standardization.
}
\value{
A numeric vector named by sample column names when available.
Samples with no usable CpG z-scores receive \code{NA}.
}
\examples{
beta <- rbind(cg00201234 = c(0.2, 0.4, 0.6),
              cg00247334 = c(0.6, 0.4, 0.2))
calculate_EMTS(beta)
}
