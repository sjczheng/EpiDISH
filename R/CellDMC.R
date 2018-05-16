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
#' A vector of phenotype. CellDMC can handle both of categorical and 
#' continuous/oderinal phenotypes. For categorical phenotypes, you must input 
#' factors to make sure you get the right results. 
#' 
#' @param frac.m
#' A matrix contains fractions of each cell-type. Each row labels a sample, 
#' with the same order of the columns in beta.m. Each column labels a 
#' cell-type. Column names, which are names of cell-types, are required. The 
#' rowSums of frac.m should be 1, and all values should be greater than 0 and 
#' less than 1.
#' 
#' @param mode
#' We provide two modes of CellDMC. One is 'basic' algorithm, and the other is 
#' 'improved' algorithm. \code{mode} can be either of 'improved' or 'basic', 
#' with default as 'improved'. For more details, pls refer to the reference.
#' 
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
#'  \code{frac.m}. For \code{basic} mode, each dataframe only contains 
#'  coefficients of DMCTs for each cell-type. For \code{improved} mode, each 
#'  dataframe contains all CpGs in input \code{beta.m}. For both modes, each 
#'  dataframe has been ranked with significant level of DMCTs of corresponding 
#'  cell-type respectively. Pls note that the order of rows in these dataframe 
#'  is different from the order of rows in \code{dmct} matrix. You can reorder 
#'  them using the rownames. All dataframes contains ranks(\code{rank}), 
#'  estimated DNAm differences(\code{Delta}), estimated T 
#'  statistics(\code{Tstat}), raw P values(\code{rawP}), and multiple 
#' hypothesis corrected P values(\code{adjP}).
#' 
#' 
#' @references 
#' Zheng SC, Beck S, Teschendorff AE. 
#' \emph{Identification of differentially methylated cell-types in 
#' Epigenome-Wide Association Studies.} In preparation (2018).
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
#' @import zoo
#' @import parallel
#' @importFrom stats model.matrix 
#' @importFrom limma lmFit eBayes topTable
#' @importFrom dplyr arrange select
#' 
#' @export
#' 
CellDMC <- function(beta.m, pheno.v, frac.m, mode = c("improved", "basic") ,
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    mc.cores = 1) {
    mode <- match.arg(mode)
    ### check input
    if (ncol(beta.m) != length(pheno.v)) 
        stop("Number of columns of beta.m should equal to length of pheno.v!")
    if (ncol(beta.m) != nrow(frac.m)) 
        stop("Number of columns of beta.m should equal to number of rows of frac.m!")
    if (length(colnames(frac.m)) != ncol(frac.m)) 
        stop("Pls assign correct name of cell-type to frac.m")
    
    ### check NA phenotype
    if (any(is.na(pheno.v))) {
      message(paste0(sum(is.na(pheno.v)), " NA phenotype found and removed."))
      retain.idx <- which(!is.na(pheno.v))
      pheno.v <- pheno.v[retain.idx]
      beta.m <- beta.m[, retain.idx]
      frac.m <- frac.m[retain.idx, ]
      if (!is.null(cov.mod)) cov.mod <- cov.mod[retain.idx, ]
    }
    
    
    ### guess factor input
    if (nlevels(factor(pheno.v)) == 2) {
      message("Binary phenotype detected. Predicted direction will be 1 - 0.")
      pheno.v <- factor(pheno.v)
    }
    if (!is.factor(pheno.v)) 
      message("pheno.v is not factor. Treating as continuous variables. Input factos for categorical phenotypes.")
    
    if (!mode %in% c("improved", "basic")) 
        stop("Input a valid mode!")
    oldw <- getOption("warn")
    options(warn = -1)
    sink("/dev/null")
    if (mode == "improved") {
      out.o <- CellDMC.improved(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, 
                                adjPMethod = adjPMethod, adjPThresh = adjPThresh, 
                                cov.mod = cov.mod, 
                                mc.cores = mc.cores)
        
    } else if (mode == "basic") {
        out.o <- CellDMC.basic(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, 
            adjPMethod = adjPMethod, adjPThresh = adjPThresh, cov.mod = cov.mod, 
            mc.cores = mc.cores)
    }
    options(warn = oldw)
    sink()
    return(out.o)
}




CellDMC.improved <- function(beta.m, pheno.v, frac.m, adjPMethod = "fdr", adjPThresh = 0.05, 
                             cov.mod = NULL, mc.cores = 1) {
    
    
    ### Fit Int models for each cell-type
    tmp.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
      design1 <- cbind(model.matrix(~pheno + pheno:frac, 
                                    data = data.frame(pheno = pheno.v, 
                                                      frac = (1 - frac.m[, j]))),
                       frac.m[, seq_len(ncol(frac.m))[-1]])
      
        if (!is.null(cov.mod)) {
          design1 <- cbind(design1, cov.mod[,-1])
        } 
        
        fit1 <- suppressWarnings(eBayes(suppressWarnings(lmFit(beta.m, design = design1))))
        m1.df <- suppressWarnings(topTable(fit1, coef = 2, number = Inf, sort.by = "none", 
            adjust.method = adjPMethod))
        return(list(m1 = m1.df))
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, ncol(frac.m)), mc.allow.recursive = TRUE, 
        mc.silent = TRUE)
    
    Delta.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$logFC))
    rawP.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$P.Value))
    adjP.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$adj.P.Val))
    T.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$t))

    
    ### pattern matching
    dmct.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
      tmpCoe.m <- rbind(Delta.m[i, ], adjP.m[i, ])
        
      sig.idx <- which(apply(tmpCoe.m, 2, function(x) {
        x[2] < adjPThresh
      }))
      ### call DMCT
      if (length(sig.idx) == 0) {
        return(c(0, rep(0, ncol(tmpCoe.m))))
      } else {
        dmct.v <- c(1, rep(0, ncol(tmpCoe.m)))
        dmct.v[sig.idx + 1] <- sign(tmpCoe.m[1, sig.idx])
        return(dmct.v)
      }

        
    }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
    
    colnames(dmct.m) <- c("DMC", colnames(frac.m))
    rownames(dmct.m) <- rownames(beta.m)
    
    ### fix overstimation
    Delta.m[Delta.m < -1] <- -1
    Delta.m[Delta.m > 1] <- 1
    
    ### coe list
    coe.ld <- lapply(seq_len(ncol(frac.m)), function(i) {
        out <- data.frame(rank = rank(rawP.m[, i], ties.method = "random"), Delta = Delta.m[, 
            i], Tstat = T.m[, i], rawP = rawP.m[, i], adjP = adjP.m[, i], cpg = rownames(beta.m))
        out <- arrange(out, rank)
        cpg <- out$cpg
        out <- as.data.frame(select(out, -cpg))
        rownames(out) <- cpg
        return(out)
    })
    names(coe.ld) <- colnames(frac.m)
    
    
    return(list(dmct = dmct.m, coe = coe.ld))
}




CellDMC.basic <- function(beta.m, pheno.v, frac.m, adjPMethod = "fdr", adjPThresh = 0.05, 
                          cov.mod = NULL, mc.cores = 1) {
    
    
    ### Fit Int models for each cell-type
    tmp.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
        design1 <- cbind(model.matrix(~pheno + pheno:frac, data = data.frame(pheno = pheno.v, 
            frac = frac.m[, j])), frac.m[, seq_len(ncol(frac.m))[-1]])
        if (!is.null(cov.mod)) {
          design1 <- cbind(design1, cov.mod[,-1])
        } 
        if (is.factor(pheno.v)) {
          design1 <- design1[, c(1, 2, 4:ncol(design1), 3)]
        }
        fit1 <- suppressWarnings(eBayes(suppressWarnings(lmFit(beta.m, design = design1))))
        m1.df <- suppressWarnings(topTable(fit1, coef = 2, number = Inf, sort.by = "none", 
            adjust.method = adjPMethod))
        m2.df <- suppressWarnings(topTable(fit1, coef = 3, number = Inf, sort.by = "none", 
            adjust.method = adjPMethod))
        return(list(m1 = m1.df, m2 = m2.df))
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, ncol(frac.m)), mc.allow.recursive = TRUE, 
        mc.silent = TRUE)
    
    
    DeltaB2.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$logFC))
    DeltaB3.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$logFC))
    Delta.m <- DeltaB2.m + DeltaB3.m
    
    rawPB2.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$P.Value))
    rawPB3.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$P.Value))
    rawP.a <- array(dim = c(nrow(rawPB2.m), ncol(rawPB2.m), 2))
    rawP.a[, , 1] <- rawPB2.m
    rawP.a[, , 2] <- rawPB3.m
    rawP.m <- apply(rawP.a, c(1, 2), min)
    
    adjPB2.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$adj.P.Val))
    adjPB3.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$adj.P.Val))
    adjP.a <- array(dim = c(nrow(adjPB2.m), ncol(adjPB2.m), 2))
    adjP.a[, , 1] <- adjPB2.m
    adjP.a[, , 2] <- adjPB3.m
    adjP.m <- apply(adjP.a, c(1, 2), min)
    
    
    
    TB2.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$t))
    TB3.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$t))
    T.a <- array(dim = c(nrow(TB2.m), ncol(TB2.m), 2))
    T.a[, , 1] <- TB2.m
    T.a[, , 2] <- TB3.m
    T.m <- apply(T.a, c(1, 2), function(x) x[which.max(abs(x))])
    
    rm(list = c("tmp.ld", "DeltaB2.m", "DeltaB3.m", "rawPB2.m", "rawPB3.m", "rawP.a", 
        "adjP.a", "TB2.m", "TB3.m", "T.a"))
    
    
    ### DMCT decision
    dmct.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
        sig.idx <- which(sapply(seq_len(ncol(frac.m)), function(j) {
            adjPB2.m[i, j] > adjPB3.m[i, j] & adjPB3.m[i, j] < adjPThresh
        }))
        
        if (length(sig.idx) == 0) {
            ### check all pattern
            if (all(adjPB2.m[i, ] < adjPThresh)) {
                return(sign(Delta.m[i, ]))
            } else {
                return(rep(0, ncol(frac.m)))
            }
            
        } else {
            dmct.v <- rep(0, ncol(frac.m))
            dmct.v[sig.idx] <- sign(Delta.m[i, sig.idx])
            return(dmct.v)
        }
        
    }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
    
    dmc.v <- apply(dmct.m, 1, function(x) ifelse(any(x != 0), 1, 0))
    dmct.m <- cbind(dmc.v, dmct.m)
    colnames(dmct.m) <- c("DMC", colnames(frac.m))
    rownames(dmct.m) <- rownames(beta.m)
    
    ### fix overstimation
    Delta.m[Delta.m < -1] <- -1
    Delta.m[Delta.m > 1] <- 1
    
    ### coe list
    coe.ld <- lapply(seq_len(ncol(frac.m)), function(i) {
        dmct.idx <- which(dmct.m[, i + 1] != 0)
        
        out <- data.frame(rank = rank(rawP.m[dmct.idx, i], ties.method = "random"), 
            Delta = Delta.m[dmct.idx, i], Tstat = T.m[dmct.idx, i], rawP = rawP.m[dmct.idx, 
                i], adjP = adjP.m[dmct.idx, i], cpg = rownames(beta.m)[dmct.idx])
        out <- arrange(out, rank)
        cpg <- out$cpg
        out <- as.data.frame(select(out, -cpg))
        rownames(out) <- cpg
        return(out)
    })
    names(coe.ld) <- colnames(frac.m)
    
    
    return(list(dmct = dmct.m, coe = coe.ld))
}

