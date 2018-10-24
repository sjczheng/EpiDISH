#' @title A method that allows the identification of differentially methylated 
#' cell-types and their directionality of change 
#' 
#' @aliases CellDMC
#' 
#' @description 
#' An outstanding challenge of Epigenome-Wide Association Studies performed in 
#' complex tissues is the identification of the specific cell-type(s) 
#' responsiblefor the observed differential methylation. CellDMC a novel 
#' statistical algorithm, which is able to identify not only differentially 
#' methylated positions, but also the specific cell-type(s) driving the 
#' differential methylation. 
#' 
#' @param beta.m
#' A beta value matrix with rows labeling the CpGs and columns labeling 
#' samples. This contains all CpGs, from which you want to find DMCTs.
#' 
#' @param pheno.v
#' A vector of phenotype. CellDMC can handle both of binary and 
#' continuous/oderinal phenotypes. \code{NA} is not allowed in
#' \code{pheno.v}.
#' 
#' @param frac.m
#' A matrix contains fractions of each cell-type. Each row labels a sample, 
#' with the same order of the columns in beta.m. Each column labels a 
#' cell-type. Column names, which are names of cell-types, are required. The 
#' rowSums of frac.m should be 1, and all values should be greater than 0 and 
#' less than 1.
#' 
#' @param adjPMethod
#' A method to adjust p values. The method can be any of method accepted by 
#' \code{\link{p.adjust}}.
#' 
#' @param adjPThresh
#' A numeric value, default as 0.05. This is used to call significant DMCTs. 
#' Adjusted p values less than this threshold will be picked.
#' 
#' @param cov.mod
#' A design matrix from \code{model.matrix}, which contains other covariates to
#' be adjusted. For example, input 
#' \code{model.matrix(~ geneder, data = pheno.df)} to adjust gender. Do not put
#'  cell-type fraction here!
#' 
#' @param sort
#' Default as \code{FALSE}. If \code{TRUE}, the data.frame in coe list will 
#' sorted based on p value of each CpG.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will be 
#' run simultaneously. The defatul is 1, which means no parallelization. 
#' 
#' @return A list with the following two items. 
#' 
#' @return dmct
#' A matrix tells wheter the input CpGs are DMCTs and DMCs. The first column 
#' gives whether a CpG is DMC or not. If the CpG is called as DMC, the value 
#' will be 1, otherwise it is 0. The following columns give DMCTs for each 
#' cell-types. If a CpG is DMCT, the value will be 1 (hypermethylated for case 
#' compared to control) or -1 (hypomethylated for case compared to control). 
#' Otherwise, the value is 0 (non-DMCT). The rows of this matrix are ordered as
#' the same as input \code{beta.m}. 
#' 
#' @return coe
#' This list contains several dataframe, which corresponds to each cel-type in
#'  \code{frac.m}. Each dataframe contains all CpGs in input \code{beta.m}. 
#'  All dataframes contain estimated DNAm differences(\code{Delta}), 
#'  standard error(\code{StdError}), estimated t statistics(\code{t}), 
#'  raw P values(\code{p}), and multiple hypothesis corrected P 
#'  values(\code{adjP}).
#' 
#' 
#' @references 
#' Zheng SC, Breeze CE, Beck S, Teschendorff AE. 
#' \emph{Identification of differentially methylated cell-types in 
#' Epigenome-Wide Association Studies.} Accepted (2018).
#' 
#' @examples 
#' data(centEpiFibIC.m)
#' data(DummyBeta.m)
#' out.l <- epidish(DummyBeta.m, centEpiFibIC.m, method = 'RPC')
#' frac.m <- out.l$estF
#' pheno.v <- factor(rep(c(0, 1), each = 5))
#' celldmc.o <- CellDMC(DummyBeta.m, pheno.v, frac.m) 
#' # Pls note this is faked beta value matrix
#' 
#' 
#' @import matrixStats
#' @import parallel
#' @importFrom stringr str_c
#' @import stats
#' @importFrom dplyr arrange select
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
        stop("Number of columns of beta.m should equal to number of rows of 
             frac.m!")
    if (length(colnames(frac.m)) != ncol(frac.m)) 
        stop("Pls assign correct name of cell-type to frac.m")
  
    
    
    ### guess factor input
    if (nlevels(factor(pheno.v)) == 2) {
        message("Binary phenotype detected. Predicted direction will be 1 - 0.")
        pheno.v <- factor(pheno.v)
    }
      
    if (!is.factor(pheno.v) & !is.character(pheno.v)) 
        message("pheno.v is not factor or character. Treating as continuous 
              variables. Input factors for categorical phenotypes.")

    
  
  ### Fit model
  design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  
  
  IntNames.v <- str_c(colnames(frac.m), "Pheno")
  colnames(design)[(1 + ncol(frac.m)):(2*ncol(frac.m))] <- IntNames.v 
  
  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
      beta.v <- beta.m[i, ]
      ### model
      Int.o <- lm(beta.v ~ ., data = data.frame(design))
    
      ### get coe
      IntCoe.m <- summary(Int.o)$coe[IntNames.v, ]
      IntCoe.v <- unlist(apply(IntCoe.m, 1, function(x) list(x)))
    
      names(IntCoe.v) <- NULL
      return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  coe.ld <- lapply(seq_len(ncol(frac.m)), function(j) {
        idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
        tmp.m <- allCoe.m[, idx]
        tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
        tmp.m[which(tmp.m[,1] > 1),1] <- 1
        tmp.m[which(tmp.m[,1] < -1),1] <- -1
        colnames(tmp.m) <- c("Delta", "StdError", "t", "p", "adjP")
        rownames(tmp.m) <- rownames(beta.m)
        return(data.frame(tmp.m))
  })
  names(coe.ld) <- colnames(frac.m)
  
  
  dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Delta")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", colnames(frac.m))
  rownames(dmct.m) <- rownames(beta.m)
  
  if(sort) {
      coe.ld <- lapply(coe.ld, function(x) {
      x$cpg <- rownames(x)
      x <- dplyr::arrange(x, "p")
      cpg <- x$cpg
      x <- as.data.frame(dplyr::select(x, -cpg))
      rownames(x) <- cpg
      return(x)
    })
  }
  
  
  return(list(dmct = dmct.m, coe = coe.ld))
}



#' @import limma
# legacy version

CellDMC.legacy <- function(beta.m, pheno.v, frac.m, 
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    sort = FALSE, mc.cores = 1) {
  ### check input
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v)) 
    stop("Number of columns of beta.m should equal to length of pheno.v!")
  if (ncol(beta.m) != nrow(frac.m)) 
    stop("Number of columns of beta.m should equal to number of rows of 
         frac.m!")
  if (length(colnames(frac.m)) != ncol(frac.m)) 
    stop("Pls assign correct name of cell-type to frac.m")
  
  
  
  ### guess factor input
  if (nlevels(factor(pheno.v)) == 2) {
    message("Binary phenotype detected. Predicted direction will be 1 - 0.")
    pheno.v <- factor(pheno.v)
  }
  
  if (!is.factor(pheno.v) & !is.character(pheno.v)) 
    message("pheno.v is not factor or character. Treating as continuous 
            variables. Input factors for categorical phenotypes.")
  
  
  coe.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
    design <- cbind(model.matrix(~type + type:frac, 
                                 data = data.frame(type = pheno.v, 
                                                   frac = (1 - frac.m[, j]))),
                    frac.m[, seq_len(ncol(frac.m))[-1]])
    if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
    fit <- suppressWarnings(eBayes(suppressWarnings(lmFit(beta.m, 
                                                          design = design))))
        m.df <- suppressWarnings(topTable(fit, coef = 2, number = Inf, 
                                      sort.by = "none",
                                       adjust.method = adjPMethod))
        coe.df <- m.df[, c("logFC", "t", "P.Value", "adj.P.Val")]
        colnames(coe.df) <- c("Delta", "t", "p", "adjP")
        rownames(coe.df) <- rownames(beta.m)
        return(coe.df)
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, ncol(frac.m)), 
     mc.allow.recursive = TRUE,
    mc.silent = TRUE)
  
    names(coe.ld) <- colnames(frac.m)
  
  
    dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
    dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
    dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Delta")[dmct.idx])
    dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
    dmct.m <- cbind(dmc.v, dmct.m)
    colnames(dmct.m) <- c("DMC", colnames(frac.m))
    rownames(dmct.m) <- rownames(beta.m)
  
    if(sort) {
        coe.ld <- lapply(coe.ld, function(x) {
        x$cpg <- rownames(x)
        x <- dplyr::arrange(x, "p")
        cpg <- x$cpg
        x <- as.data.frame(dplyr::select(x, -cpg))
        rownames(x) <- cpg
        return(x)
    })
  }
  
  
    return(list(dmct = dmct.m, coe = coe.ld))
}

