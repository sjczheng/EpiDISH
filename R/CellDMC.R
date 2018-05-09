#' @title A method that allows the identification of differentially methylated 
#' cell-types and their directionality of change 
#' 
#' @aliases CellDMC
#' 
#' @description 
#' An outstanding challenge of Epigenome-Wide Association Studies performed in 
#' complex tissues is the identification of the specific cell-type(s) responsible
#'  for the observed differential methylation. CellDMC a novel statistical 
#'  algorithm, which is able to identify not only differentially methylated positions,
#'   but also the specific cell-type(s) driving the differential methylation. 
#' 
#' @param beta.m
#' A beta value matrix with rows labeling the CpGs and columns labeling samples. 
#' This contains all CpGs, from which you want to find DMCTs.
#' 
#' @param pheno.v
#' A vector of phenotype. For binary phenotypes, each item has to be 0 or 1, 
#' with 1 labeling case and 0 labeling control.
#' 
#' @param frac.m
#' A matrix contains fractions of each cell-type. Each row labels a sample, with 
#' the same order of the columns in beta.m. Each column labels a cell-type. 
#' Column names, which are names of cell-types, are required. The rowSums of 
#' frac.m should be 1, and all values should be greater than 0 and less than 1.
#' 
#' @param mode
#' We provide two modes of CellDMC. One is 'basic' algorithm, and the other is '
#' improved' algorithm. \code{mode} can be either of 'improved' or 'basic', with 
#' default as 'improved'. For more details, pls refer to the reference.
#' 
#' @param pheno.class
#' A string tell CellDMC the class of phenotypes, which can be either of 
#' \code{bi}(for binary phenotypes) or \code{continuous}(for continuous phenotypes).
#' 
#' @param adjPMethod
#' A method to adjust p values. The method can be any of method accepted by 
#' \code{\link{p.adjust}}.
#' 
#' @param adjPThresh
#' A numeric value, default as 0.05. This is used to call significant DMCTs. 
#' Adjusted p values less than this threshold will be picked.
#' 
#' @param DiffThresh
#' A DNAm diff threshold. The default is 0.1. For each cell-type, CpGs with 
#' absolute DNAm change greater than this threshold will be treated as DMCTs. 
#' Pls note that this threshold is only used in 'improved' mode for binary pheno.
#' 
#' @param mc.cores
#' The number of cores to use, i.e. at most how many child processes will be run 
#' simultaneously. The defatul is 1, which means no parallelization. 
#' 
#' @return A list with the following two items. 
#' 
#' @return dmct
#' A matrix tells wheter the input CpGs are DMCTs and DMCs. The first column gives 
#' whether a CpG is DMC or not. If the CpG is called as DMC, the value will be 1,
#' otherwise it is 0. The following columns give DMCTs for each cell-types. If a 
#' CpG is DMCT, the value will be 1 (hypermethylated for case compared to control) 
#' or -1 (hypomethylated for case compared to control). Otherwise, the value is 0 
#' (non-DMCT). The rows of this matrix are ordered as the same as input 
#' \code{beta.m}. 
#' 
#' @return coe
#' This list contains several dataframe, which corresponds to each cel-type in \code{frac.m}.
#' For \code{basic} mode, each dataframe only contains coefficients of DMCTs for 
#' each cell-type. For \code{improved} mode, each dataframe contains all CpGs in 
#' input \code{beta.m}. For both modes, each dataframe has been ranked with 
#' significant level of DMCTs of corresponding cell-type respectively. Pls note 
#' that the order of rows in these dataframe is different from the order of rows 
#' in \code{dmct} matrix. You can reorder them using the rownames. All dataframes 
#' contains ranks(\code{rank}), estimated DNAm differences(\code{Delta}), 
#' estimated T statistics(\code{Tstat}), raw P values(\code{rawP}), and multiple 
#' hypothesis corrected P values(\code{adjP}).
#' 
#' 
#' @references 
#' Zheng SC, Beck S, Teschendorff AE. 
#' \emph{Identification of differentially methylated cell-types in Epigenome-Wide 
#' Association Studies.} In preparation (2018).
#' 
#' @examples 
#' data(centEpiFibIC.m)
#' data(DummyBeta.m)
#' out.l <- epidish(DummyBeta.m, centEpiFibIC.m, method = 'RPC')
#' frac.m <- out.l$estF
#' pheno.v <- rep(c(0, 1), each = 5)
#' celldmc.o <- CellDMC(DummyBeta.m, pheno.v, frac.m) # Pls note this is faked beta value matrix
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
CellDMC <- function(beta.m, pheno.v, frac.m, mode = c("improved", "basic"), pheno.class = c("bi", "continuous") ,adjPMethod = "fdr", 
    adjPThresh = 0.05, DiffThresh = 0.1, mc.cores = 1) {
    mode <- match.arg(mode)
    pheno.class <- match.arg(pheno.class)
    ### check input
    if (ncol(beta.m) != length(pheno.v)) 
        stop("Number of columns of beta.m should equal to length of pheno.v!")
    if (ncol(beta.m) != nrow(frac.m)) 
        stop("Number of columns of beta.m should equal to number of rows of frac.m!")
    if (length(colnames(frac.m)) != ncol(frac.m)) 
        stop("Pls assign correct name of cell-type to frac.m")
    if (pheno.class == "bi") {
      if (sum(pheno.v %in% c(0, 1)) != length(pheno.v)) 
        stop("Pls code case as 1 and control as 0 in pheno.v. CellDMC only accepts binary phenotypes or continuous phenotypes for now.")
    }

    
    pheno.v <- factor(as.character(pheno.v), levels = as.character(c(0, 1)))
    
    if (!mode %in% c("improved", "basic")) 
        stop("Input a valid mode!")
    if (!pheno.class %in% c("bi", "continuous")) 
      stop("Input a valid pheno.class!")
    oldw <- getOption("warn")
    options(warn = -1)
    sink("/dev/null")
    if (mode == "improved") {
      if (pheno.class == "bi") {
        out.o <- CellDMC.improved(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, 
                                  adjPMethod = adjPMethod, adjPThresh = adjPThresh, DiffThresh = DiffThresh, 
                                  mc.cores = mc.cores)
      } else if (pheno.class == "continuous"){
        out.o <- CellDMC.continuous(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, 
                                    adjPMethod = adjPMethod, adjPThresh = adjPThresh, 
                                    mc.cores = mc.cores)
      }
        
    } else if (mode == "basic") {
        out.o <- CellDMC.basic(beta.m = beta.m, pheno.v = pheno.v, frac.m = frac.m, 
            adjPMethod = adjPMethod, adjPThresh = adjPThresh, mc.cores = mc.cores)
    }
    options(warn = oldw)
    sink()
    return(out.o)
}




CellDMC.improved <- function(beta.m, pheno.v, frac.m, adjPMethod = "fdr", adjPThresh = 0.05, 
    DiffThresh = 0.1, mc.cores = 1) {
    
    ### Fit NoInt model
    design0 <- cbind(model.matrix(~type, data = data.frame(type = pheno.v)), frac.m[, 
        seq_len(ncol(frac.m))[-1]])
    fit0 <- eBayes(lmFit(beta.m, design = design0))
    m0.adjP <- topTable(fit0, coef = 2, number = Inf, sort.by = "none", adjust.method = adjPMethod)$adj.P.Val
    
    ### Fit Int models for each cell-type
    tmp.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
        design1 <- cbind(model.matrix(~type + type:frac, data = data.frame(type = pheno.v, 
            frac = (1 - frac.m[, j]))), frac.m[, seq_len(ncol(frac.m))[-1]])
        fit1 <- suppressWarnings(eBayes(suppressWarnings(lmFit(beta.m, design = design1))))
        m1.df <- suppressWarnings(topTable(fit1, coef = 2, number = Inf, sort.by = "none", 
            adjust.method = adjPMethod))
        design2 <- cbind(model.matrix(~type + type:frac, data = data.frame(type = pheno.v, 
            frac = frac.m[, j])), frac.m[, seq_len(ncol(frac.m))[-1]])
        fit2 <- suppressWarnings(eBayes(suppressWarnings(lmFit(beta.m, design = design2))))
        m2.df <- suppressWarnings(topTable(fit2, coef = 2, number = Inf, sort.by = "none", 
            adjust.method = adjPMethod))
        return(list(m1 = m1.df, m2 = m2.df))
    }, mc.preschedule = TRUE, mc.cores = min(mc.cores, ncol(frac.m)), mc.allow.recursive = TRUE, 
        mc.silent = TRUE)
    
    Delta.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$logFC))
    rawP.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$P.Value))
    adjP.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$adj.P.Val))
    T.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m1$t))
    tmpDelta.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$logFC))
    tmpAdjP.m <- do.call(cbind, lapply(tmp.ld, function(x) x$m2$adj.P.Val))
    tmpSign.m <- do.call(cbind, lapply(tmp.ld, function(x) sign(x$m1$logFC) * sign(x$m2$logFC)))
    rm(tmp.ld)
    
    is.uniAll <- function(tmpCoe.m, i) {
        all(tmpCoe.m[5, ] == 1) & (nlevels(factor(sign(tmpCoe.m[1, ]))) == 1) & all(apply(tmpCoe.m, 
            2, function(x) all(abs(x[c(1, 3)]) > DiffThresh))) & m0.adjP[i] < adjPThresh
    }
    is.Bi <- function(tmpCoe.m) {
        any(apply(tmpCoe.m, 2, function(x) {
            abs(x[1]) > DiffThresh & abs(x[3]) > DiffThresh & x[2] < adjPThresh & 
                x[4] < adjPThresh & x[5] == -1
        })) & nlevels(factor(sign(tmpCoe.m[1, ]))) == 2
        
    }
    
    ### pattern matching
    dmct.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
        tmpCoe.m <- rbind(Delta.m[i, ], adjP.m[i, ], tmpDelta.m[i, ], tmpAdjP.m[i, 
            ], tmpSign.m[i, ])
        
        ### check uni-all pattern
        if (DiffThresh > 0) {
            if (is.uniAll(tmpCoe.m, i)) 
                return(c(1, rep(sign(tmpCoe.m[1, 1]), ncol(tmpCoe.m))))
        }
        
        sig.idx <- which(apply(tmpCoe.m, 2, function(x) {
            abs(x[1]) > DiffThresh & x[2] < adjPThresh
        }))
        
        ### check bi-directional pattern
        if (nlevels(factor(sign(tmpCoe.m[1, sig.idx]))) == 1) {
            if (is.Bi(tmpCoe.m)) {
                tmp.m <- cbind(tmpCoe.m[1, ], tmpCoe.m[2, ], seq_len(ncol(frac.m)))
                tmp.m <- tmp.m[which(sign(tmp.m[, 1]) != sign(tmp.m[sig.idx[1], 1])), 
                  , drop = FALSE]
                sig.idx <- c(sig.idx, tmp.m[which.min(tmp.m[, 2]), 3])
                dmc.v <- c(1, rep(0, ncol(tmpCoe.m)))
                dmc.v[sig.idx + 1] <- sign(tmpCoe.m[1, sig.idx])
                return(dmc.v)
            }
        }
        
        ### all other
        if (length(sig.idx) == 0) {
            return(c(0, rep(0, ncol(tmpCoe.m))))
        } else {
            dmc.v <- c(1, rep(0, ncol(tmpCoe.m)))
            dmc.v[sig.idx + 1] <- sign(tmpCoe.m[1, sig.idx])
            return(dmc.v)
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


CellDMC.continuous <- function(beta.m, pheno.v, frac.m, adjPMethod = "fdr", adjPThresh = 0.05, 
                              mc.cores = 1) {
  
  ### Fit Int models for each cell-type
  tmp.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
    design1 <- cbind(model.matrix(~pheno + pheno:frac, data = data.frame(pheno = pheno.v, 
                                                                       frac = (1 - frac.m[, j]))), frac.m[, seq_len(ncol(frac.m))[-1]])
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
  
    if (length(sig.idx) == 0) {
      return(c(0, rep(0, ncol(tmpCoe.m))))
    } else {
      dmc.v <- c(1, rep(0, ncol(tmpCoe.m)))
      dmc.v[sig.idx + 1] <- sign(tmpCoe.m[1, sig.idx])
      return(dmc.v)
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
    mc.cores = 1) {
    
    
    ### Fit Int models for each cell-type
    tmp.ld <- mclapply(seq_len(ncol(frac.m)), function(j) {
        design1 <- cbind(model.matrix(~type + type:frac, data = data.frame(type = pheno.v, 
            frac = frac.m[, j])), frac.m[, seq_len(ncol(frac.m))[-1]])
        design1 <- design1[, c(1, 2, 4:ncol(design1), 3)]
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
            dmc.v <- rep(0, ncol(frac.m))
            dmc.v[sig.idx] <- sign(Delta.m[i, sig.idx])
            return(dmc.v)
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

