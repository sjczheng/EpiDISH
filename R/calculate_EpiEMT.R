# Input: numeric beta matrix with CpG row names and samples in columns.
# Output: numeric EMTS vector named by sample column names, when supplied.
# Each CpG is standardized across the supplied cohort before signed averaging.
# Missing CpGs and non-finite values are skipped. CpGs with fewer than two
# finite observations or zero variance are excluded; no usable CpGs gives NA.
# Usage: scores <- calculate_EMTS(beta)

calculate_EMTS <- function(beta) {
  # Fixed age-effect directions from score_cpgs.csv (33 CpGs, 21 EMT genes).
  direction <- c(
    cg00201234 = 1, cg00247334 = -1, cg02271895 = 1, cg02327997 = 1,
    cg02935338 = 1, cg04081201 = 1, cg05237641 = 1, cg06266514 = 1,
    cg06375949 = 1, cg07035659 = 1, cg07249459 = 1, cg08009622 = 1,
    cg09255732 = -1, cg09949775 = 1, cg10266211 = 1, cg11668399 = 1,
    cg12681727 = 1, cg14621053 = 1, cg15574301 = 1, cg16670307 = 1,
    cg16695483 = 1, cg17285225 = 1, cg18390610 = 1, cg19214331 = 1,
    cg19485202 = 1, cg19907305 = 1, cg20640433 = 1, cg22670329 = 1,
    cg22881914 = 1, cg25144207 = 1, cg26669793 = 1, cg26948907 = 1,
    cg26949694 = 1
  )

  if (!is.matrix(beta) || !is.numeric(beta) || is.null(rownames(beta))) {
    stop("beta must be a numeric matrix with CpG row names and samples in columns.")
  }
  cpgs <- intersect(names(direction), rownames(beta))
  x <- beta[cpgs, , drop = FALSE]
  scores <- setNames(rep(NA_real_, ncol(x)), colnames(x))
  if (nrow(x) == 0L || ncol(x) < 2L) return(scores)
  x[!is.finite(x)] <- NA_real_

  # Original definition: standardize each CpG within the supplied cohort.
  cpg_mean <- rowMeans(x, na.rm = TRUE)
  cpg_sd <- apply(x, 1L, sd, na.rm = TRUE)
  usable <- is.finite(cpg_sd) & cpg_sd > 0
  if (!any(usable)) return(scores)
  z <- sweep(x[usable, , drop = FALSE], 1L, cpg_mean[usable], "-")
  z <- sweep(z, 1L, cpg_sd[usable], "/")
  z <- sweep(z, 1L, direction[cpgs[usable]], "*")

  # Average each sample's available direction-adjusted z-scores.
  z[!is.finite(z)] <- NA_real_
  scores <- colMeans(z, na.rm = TRUE)
  scores[!is.finite(scores)] <- NA_real_
  return(scores)
}
