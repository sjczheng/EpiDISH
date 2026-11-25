# Input: EpiDISH::epidish() output list (containing $estF) or a numeric cell-type
#        fraction matrix/data.frame with samples in rows and cell types in columns
#        (must include "Fibroblast" and "Epithelial").
# Output: numeric mFDI vector named by sample row names, when supplied.
# Formula: log2((Fibroblast + 0.01) / (Epithelial + 0.01))
# Usage:
#   est <- EpiDISH::epidish(beta.m = beta, ref.m = ref, method = "RPC")
#   mfdi <- calculate_mFDI(est$estF)  # or calculate_mFDI(est)

calculate_mFDI <- function(estF, epsilon = 0.01) {
  if (is.list(estF) && !is.data.frame(estF) && "estF" %in% names(estF)) {
    estF <- estF$estF
  }
  if (!(is.matrix(estF) || is.data.frame(estF)) || is.null(colnames(estF))) {
    stop("estF must be an EpiDISH::epidish() result or a matrix/data.frame with cell-type column names.")
  }

  required_cols <- c("Fibroblast", "Epithelial")
  missing_cols <- setdiff(required_cols, colnames(estF))
  if (length(missing_cols) > 0L) {
    stop("Missing required cell-type column(s): ", paste(missing_cols, collapse = ", "))
  }

  Fibroblast <- as.numeric(estF[, "Fibroblast"])
  epithelial <- as.numeric(estF[, "Epithelial"])

  mfdi <- rep(NA_real_, nrow(estF))
  valid <- is.finite(Fibroblast) & is.finite(epithelial) &
    (Fibroblast + epsilon > 0) & (epithelial + epsilon > 0)
  mfdi[valid] <- log2((Fibroblast[valid] + epsilon) / (epithelial[valid] + epsilon))

  if (!is.null(rownames(estF))) {
    names(mfdi) <- rownames(estF)
  }
  return(mfdi)
}
