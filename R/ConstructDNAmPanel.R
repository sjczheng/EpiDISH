#' @title Construct DNA Methylation reference Panel
#'
#' @description Constructs a cell-type-specific DNA methylation panel by selecting
#'              top-ranked hypo-methylated probes based on a "Gap Specificity Score" metric,
#'              which prioritizes markers with large effect sizes.
#'
#' @param dnam.matrix A DNA methylation beta value matrix with no NA values.
#'                    Rows should be Probe IDs and columns should be samples.
#' @param cell.types.info A data frame that must contain a 'CellType' column,
#'                        specifying the cell type for each sample in `dnam.matrix`.
#' @param panel.cell.types A character vector specifying for which cell types a
#'                         panel should be constructed. These must be present in
#'                         `cell.types.info$CellType`.
#' @param n.probes The number of top-ranked hypo-methylated probes to select for
#'                 each cell type. Default is 100.
#' @param fdr.threshold The FDR (False Discovery Rate) threshold for filtering
#'                      significantly differentially methylated probes. Default is 0.05.
#' @param p.adjust.method The method for multiple testing correction. Must be one of
#'                        the methods supported by `stats::p.adjust`. Default is "BH".
#' @param summary.method The method to use for summarizing methylation values when
#'                       creating the reference panel. Can be "mean" or "median".
#'                       Default is "mean".
#' @param equal.variance A logical value (`TRUE` or `FALSE`) passed to `genefilter::fastT`, 
#'                       indicating whether to assume equal variance between groups. Default is FALSE.
#'
#' @return A list containing two components:
#'         \item{DetailedResults}{A list where each element is a data frame for a
#'                                cell type in `panel.cell.types`. The data frame
#'                                contains detailed information for the selected probes,
#'                                such as GapScore, AbsEffectSize, Zstatistic, and Pvalue.}
#'         \item{ReferencePanel}{A matrix with probes as rows and cell types as columns.
#'                               The values are the mean or median methylation beta values
#'                               for each probe in the corresponding cell type. This can serve
#'                               as a reference matrix for deconvolution.}
#'
#' @details
#' The function's logic is as follows:
#' 1. For each target cell type, it is compared against all other cell types.
#' 2. A t-test is used to find differentially methylated probes, and p-values are
#'    adjusted to FDR using the specified `p.adjust.method`.
#' 3. Probes are filtered based on the `fdr.threshold`.
#' 4. Only probes that are hypo-methylated in the target cell type (negative t-statistic)
#'    are retained.
#' 5. A "Gap specificity Score" is calculated for these probes (defined as the minimum beta value
#'    in other cell types minus the maximum beta value in the target cell type).
#' 6. Probes are ranked by the Gap Score, and the top `n.probes` are selected.
#' 7. Finally, the probes selected for all target cell types are combined to form the
#'    final panel, and the reference matrix is computed using the specified `summary.method`.
#'
#' @importFrom genefilter fastT
#' @importFrom stats p.adjust pt p.adjust.methods
#' @importFrom matrixStats rowMaxs rowMins rowMeans2 rowMedians
#' @importFrom utils head
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have referenceBetas_matrix_noNA and referenceCovars_df
#' # library(genefilter)
#' # library(stats)
#' # library(matrixStats)
#' 
#' # Define target cell types for the panel
#' target_cell_types <- c("CD4T", "Bcell", "Mono")
#' 
#' # Construct the Panel
#' panel_results <- ConstructDNAmPanel(
#'   dnam.matrix = referenceBetas_matrix_noNA,
#'   cell.types.info = referenceCovars_df,
#'   panel.cell.types = target_cell_types,
#'   n.probes = 50,
#'   fdr.threshold = 0.05,
#'   p.adjust.method = "BH",
#'   summary.method = "median" 
#' )
#' 
#' # View the results
#' # View the reference panel matrix
#' head(panel_results$ReferencePanel)
#' 
#' # View detailed probe info for B-cells
#' head(panel_results$DetailedResults$Bcell)
#' }
ConstructDNAmPanel <- function(dnam.matrix, 
                               cell.types.info, 
                               panel.cell.types, 
                               n.probes = 100, 
                               fdr.threshold = 0.05, 
                               p.adjust.method = "BH",
                               summary.method = "mean",
                               equal.variance = FALSE) {
  
  # 1. Check for required packages
  if (!requireNamespace("genefilter", quietly = TRUE)) {
    stop("Package 'genefilter' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("Package 'stats' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    warning("Package 'matrixStats' is recommended for performance.", call. = FALSE)
    rowMeans2 <- function(x, ...) apply(x, 1, mean, ...)
    rowMaxs <- function(x, ...) apply(x, 1, max, ...)
    rowMins <- function(x, ...) apply(x, 1, min, ...)
  } else {
    rowMeans2 <- matrixStats::rowMeans2
    rowMaxs <- matrixStats::rowMaxs
    rowMins <- matrixStats::rowMins
    rowMedians <- matrixStats::rowMedians
  }
  
  # 2. Validate inputs
  if (!is.matrix(dnam.matrix)) dnam.matrix <- as.matrix(dnam.matrix)
  if (anyNA(dnam.matrix)) stop("'dnam.matrix' must be NA-free.")
  if (!"CellType" %in% colnames(cell.types.info)) stop("'cell.types.info' must have a 'CellType' column.")
  if (nrow(cell.types.info) != ncol(dnam.matrix)) stop("Row count of 'cell.types.info' must match column count of 'dnam.matrix'.")
  if (is.null(rownames(dnam.matrix))) stop("'dnam.matrix' must have probe IDs as rownames.")
  if (!all(panel.cell.types %in% cell.types.info$CellType)) {
    stop("Some 'panel.cell.types' not found in 'cell.types.info$CellType'.")
  }
  if (!summary.method %in% c("mean", "median")) {
    stop("'summary.method' must be either 'mean' or 'median'.")
  }
  if (summary.method == "median" && !requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package 'matrixStats' is required to use summary.method = 'median'.")
  }
  # --- New Validation Rules ---
  if (!is.logical(equal.variance) || length(equal.variance) != 1) {
    stop("'equal.variance' must be a single logical value (TRUE or FALSE).")
  }
  if (!p.adjust.method %in% stats::p.adjust.methods) {
    stop(paste("'p.adjust.method' is not valid. Choose one of:", 
               paste(stats::p.adjust.methods, collapse = ", ")))
  }
  # --------------------------
  
  message("Input validation complete.")
  
  # 3. Preprocess data
  tIndexes_all <- split(seq_len(ncol(dnam.matrix)), cell.types.info$CellType)
  
  cell_types_present_in_data <- intersect(panel.cell.types, names(tIndexes_all))
  if (length(cell_types_present_in_data) < 2) {
    stop("Fewer than two cell types are available to build the panel.")
  }
  
  message("Data preprocessing complete.")
  
  # 4. Core logic: select probes for each cell type
  detailed_results_list <- vector("list", length(cell_types_present_in_data))
  names(detailed_results_list) <- cell_types_present_in_data
  all_selected_probes <- character(0)
  
  for (target_ct_name in cell_types_present_in_data) {
    message("Processing cell type: ", target_ct_name)
    
    x1 <- tIndexes_all[[target_ct_name]]
    other_cts_names <- setdiff(names(tIndexes_all), target_ct_name)
    x2_all_others <- unlist(tIndexes_all[other_cts_names])
    
    if (length(x1) < 1 || length(x2_all_others) < 1) {
      warning(paste("Skipping", target_ct_name, "- zero samples in target or other groups."))
      next
    }
    
    df <- length(x1) + length(x2_all_others) - 2
    if (df <= 0) {
      warning(paste("Skipping", target_ct_name, "- non-positive degrees of freedom."))
      next
    }
    
    tstat_result <- tryCatch(genefilter::fastT(dnam.matrix, x1, x2_all_others, var.equal = equal.variance),
                             error = function(e) { warning(paste("t-test failed for", target_ct_name, ":", e$message)); NULL })
    
    if (is.null(tstat_result)) next
    
    p_values <- 2 * stats::pt(abs(tstat_result$z), df = df, lower.tail = FALSE)
    fdr <- stats::p.adjust(p_values, method = p.adjust.method)
    
    sig_probes_idx <- which(fdr < fdr.threshold)
    if (length(sig_probes_idx) == 0) {
      message("No significant probes for ", target_ct_name, " at FDR < ", fdr.threshold)
      next
    }
    
    p_sig <- dnam.matrix[sig_probes_idx, , drop = FALSE]
    tstat_sig_z <- tstat_result$z[sig_probes_idx]
    
    is_hypo <- sign(tstat_sig_z) < 0
    hypo_probes_idx_in_sig <- which(is_hypo)
    
    if (length(hypo_probes_idx_in_sig) == 0) {
      message("No hypo-methylated probes found for ", target_ct_name)
      next
    }
    
    max_target_vals <- rowMaxs(p_sig[hypo_probes_idx_in_sig, x1, drop = FALSE])
    min_other_samples <- rowMins(p_sig[hypo_probes_idx_in_sig, x2_all_others, drop = FALSE])
    gap_score_hypo <- min_other_samples - max_target_vals
    
    hypo_probes_in_sig <- rownames(p_sig)[hypo_probes_idx_in_sig]
    sorted_hypo_probes <- hypo_probes_in_sig[order(gap_score_hypo, decreasing = TRUE)]
    
    num_to_select <- min(n.probes, length(sorted_hypo_probes))
    final_selected_probes <- head(sorted_hypo_probes, num_to_select)
    
    if (length(final_selected_probes) == 0) next
    
    all_selected_probes <- c(all_selected_probes, final_selected_probes)
    
    selected_indices_in_sig <- match(final_selected_probes, rownames(p_sig))
    mean_target <- rowMeans2(p_sig[selected_indices_in_sig, x1, drop = FALSE])
    mean_all_others <- rowMeans2(p_sig[selected_indices_in_sig, x2_all_others, drop = FALSE])
    
    df_ct <- data.frame(
      ProbeID = final_selected_probes,
      GapScore = gap_score_hypo[match(final_selected_probes, hypo_probes_in_sig)],
      AbsEffectSize = abs(mean_target - mean_all_others),
      Zstatistic = tstat_sig_z[selected_indices_in_sig],
      Pvalue = p_values[sig_probes_idx][selected_indices_in_sig],
      row.names = final_selected_probes
    )
    detailed_results_list[[target_ct_name]] <- df_ct
  }
  
  # 5. Finalize results
  detailed_results_list <- detailed_results_list[!sapply(detailed_results_list, is.null)]
  candidatePanelProbes <- unique(all_selected_probes)
  
  if (length(candidatePanelProbes) == 0) {
    warning("No probes were selected for the final panel.")
    return(list(DetailedResults = detailed_results_list, ReferencePanel = matrix(0,0,0)))
  }
  
  # 6. Create the reference panel matrix
  summary.func <- if (summary.method == "mean") rowMeans2 else rowMedians
  
  p_candidate <- dnam.matrix[candidatePanelProbes, , drop = FALSE]
  processed_cell_types <- names(detailed_results_list)
  
  ReferencePanel <- matrix(NA, nrow = length(candidatePanelProbes), ncol = length(processed_cell_types))
  rownames(ReferencePanel) <- candidatePanelProbes
  colnames(ReferencePanel) <- processed_cell_types
  
  for (ct_name in processed_cell_types) {
    ind <- tIndexes_all[[ct_name]]
    if (length(ind) > 0) {
      ReferencePanel[, ct_name] <- summary.func(p_candidate[, ind, drop = FALSE])
    }
  }
  
  message("Panel construction complete. Returning results.")
  
  return(list(DetailedResults = detailed_results_list, 
              ReferencePanel = ReferencePanel))
}
