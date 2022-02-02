#' @title A function that allows the identification of differentially methylated 
#' cell-types in in Epigenome-Wide Association Studies(EWAS)
#' 
#' @aliases CellDMC
#' 
#' @description 
#' An outstanding challenge of Epigenome-Wide Association Studies performed in 
#' complex tissues is the identification of the specific cell-type(s) 
#' responsible for the observed differential methylation. CellDMC is a novel 
#' statistical algorithm, which is able to identify not only differentially 
#' methylated positions, but also the specific cell-type(s) driving the 
#' methylation change. 
#' 
#' @param beta.m
#' A beta value matrix with rows labeling the CpGs and columns labeling 
#' samples. 
#' 
#' @param pheno.v
#' A vector of phenotype. CellDMC can handle both of binary and 
#' continuous/oderinal phenotypes. \code{NA} is not allowed in
#' \code{pheno.v}.
#' 
#' @param frac.m
#' A matrix contains fractions of each cell-type. Each row labels a sample, 
#' with the same order of the columns in beta.m. Each column labels a 
#' cell-type. Column names, which are the names of cell-types, are required. The 
#' rowSums of frac.m should be 1 or close to 1.
#' 
#' @param adjPMethod
#' The method used to adjust p values. The method can be any of method
#'  accepted by \code{\link{p.adjust}}.
#' 
#' @param adjPThresh
#' A numeric value, default as 0.05. This is used to call DMCTs. 
#' For each cell-type respectively, the CpG with the adjusted p values less than
#'  this threshold will be reported as DMCTs (-1 or 1) in the 'dmct' matrix in 
#'  the returned list.
#' 
#' @param cov.mod
#' A design matrix from \code{model.matrix}, which contains other covariates to
#' be adjusted. For example, input 
#' \code{model.matrix(~ geneder, data = pheno.df)} to adjust gender. Do not put
#'  cell-type fraction here!
#' 
#' @param sort
#' Default as \code{FALSE}. If \code{TRUE}, the data.frame in coe list will be
#' sorted based on p value of each CpG. The order of rows in 'dmct' will not 
#' change since the orders of each cell-type are different.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many threads will 
#' run simultaneously. The defatul is 1, which means no parallelization. 
#' 
#' @return A list with the following two items. 
#' 
#' @return dmct
#' A matrix gives wheter the input CpGs are DMCTs and DMCs. The first column 
#' tells whether a CpG is a DMC or not. If the CpG is called as DMC, the value 
#' will be 1, otherwise it is 0. The following columns give DMCTs for each 
#' cell-type. If a CpG is a DMCT, the value will be 1 (hypermethylated for case 
#' compared to control) or -1 (hypomethylated for case compared to control). 
#' Otherwise, the value is 0 (non-DMCT). The rows of this matrix are ordered as
#' the same as that of the input \code{beta.m}. 
#' 
#' @return coe
#' This list contains several dataframes, which correspond to each cell-type in
#'  \code{frac.m}. Each dataframe contains all CpGs in input \code{beta.m}. 
#'  All dataframes contain estimated DNAm changes (\code{Estimate}), 
#'  standard error (\code{SE}), estimated t statistics (\code{t}), 
#'  raw P values (\code{p}), and multiple hypothesis corrected P 
#'  values (\code{adjP}).
#' 
#' 
#' @references 
#' Zheng SC, Breeze CE, Beck S, Teschendorff AE. 
#' \emph{Identification of differentially methylated cell-types in 
#' Epigenome-Wide Association Studies.} 
#' Nat Methods (2018) 15: 1059-1066
#' doi:\href{https://doi.org/10.1038/s41592-018-0213-x}{10.1038/s41592-018-0213-x}.
#' 
#' @examples 
#' data(centEpiFibIC.m)
#' data(DummyBeta.m)
#' out.l <- epidish(DummyBeta.m, centEpiFibIC.m, method = 'RPC')
#' frac.m <- out.l$estF
#' pheno.v <- rep(c(0, 1), each = 5)
#' celldmc.o <- CellDMC(DummyBeta.m, pheno.v, frac.m) 
#' # Pls note this is a faked beta value matrix.
#' 
#' 
#' @import matrixStats
#' @import parallel
#' @importFrom stringr str_c
#' @importFrom Matrix rankMatrix
#' @import stats
#' 
#' @export
#' 
CellDMC <- function(beta.m, pheno.v, frac.m, 
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    sort = FALSE, mc.cores = 1) {
    ### check input
    if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
    if (ncol(beta.m) != length(pheno.v)) 
        stop("Number of columns of beta.m should equal to length of pheno.v!")
    if (ncol(beta.m) != nrow(frac.m)) 
        stop("Number of columns of beta.m should equal to number of rows of frac.m!")
    if (length(colnames(frac.m)) != ncol(frac.m)) 
        stop("Pls assign correct names of cell-type to frac.m")
  
    ### check whether input is beta value matrix
    is.beta <- ((min(beta.m) >= 0) & (max(beta.m) <= 1))
    
    ### guess factor input
    if (nlevels(factor(pheno.v)) == 2) {
        message("Binary phenotype detected. Predicted change will be 1 - 0.")
        pheno.v <- factor(pheno.v)
    }
    if (!is.factor(pheno.v) & !is.character(pheno.v)) 
        message("pheno.v is not factor or character. Treating as continuous variables.")

  ### Fit model
  design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
  if (Matrix::rankMatrix(design) < ncol(design)) {
    stop("The design matrix is not full ranked.\nThis means that you coundn't make inference for all cell-types in your fraction matrix.
         This is usally casued by fractions of a cell-type of one pheno type are all 0 or some fractions in one pheno type are paralle to that of another cell-type.
         You might use which(colSums(model.matrix(~ frac.m + pheno.v:frac.m)[, -1]) == 0) to find the cell type.")
  }
  
  
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  IntNames.v <- str_c(colnames(frac.m), "Pheno")
  colnames(design)[(1 + ncol(frac.m)):(2*ncol(frac.m))] <- IntNames.v 
  
  ### fit linear model for each CpG
  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
      beta.v <- beta.m[i, ]
      ### model
      Int.o <- lm(beta.v ~ .-1, data = data.frame(design))
    
      ### get coe
      IntCoe.m <- summary(Int.o)$coe[IntNames.v, ]
      IntCoe.v <- unlist(apply(IntCoe.m, 1, function(x) list(x)))
    
      names(IntCoe.v) <- NULL
      return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  ### extract coefficients for each cell-type
  coe.ld <- lapply(seq_len(ncol(frac.m)), function(j) {
        idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
        tmp.m <- allCoe.m[, idx]
        tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
        if (is.beta) { 
            tmp.m[which(tmp.m[,1] > 1),1] <- 1
            tmp.m[which(tmp.m[,1] < -1),1] <- -1
        }  ### if input is a beta values matrix, bound the estimated changes

        colnames(tmp.m) <- c("Estimate", "SE", "t", "p", "adjP")
        rownames(tmp.m) <- rownames(beta.m)
        return(data.frame(tmp.m))
  })
  names(coe.ld) <- colnames(frac.m)
  
  ### get dmct matrix
  dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", colnames(frac.m))
  rownames(dmct.m) <- rownames(beta.m)
  
  if(sort) coe.ld <- lapply(coe.ld, function(x) x[order(x$p),] )
  
  return(list(dmct = dmct.m, coe = coe.ld))
}

