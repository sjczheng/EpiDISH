#' @title 
#' Epigenetic Dissection of Intra-Sample-Heterogeneity
#' 
#' @aliases epidish
#'  
#' @description 
#' A reference-based function to infer the fractions of a priori known cell 
#' subtypes present in a sample representing a mixture of such cell-types. 
#' Inference proceeds via one of 3 methods (Robust Partial Correlations-RPC, 
#' Cibersort-CBS, Constrained Projection-CP), as determined by the user.
#' 
#' 
#' @param beta.m
#' A data matrix with rows labeling the molecular features (should use same ID 
#' as in ref.m) and columns labeling samples (e.g. primary tumour specimens). 
#' Missing value is not allowed and all values should be positive or zero. 
#' In the case of DNA methylation, these are beta-values.
#' 
#' @param ref.m
#' A matrix of reference 'centroids', i.e. representative molecular profiles, 
#' for a number of cell subtypes. rows label molecular features (e.g. CpGs,...) 
#' and columns label the cell-type. IDs need to be provided as rownames and 
#' colnames, respectively. Missing value is not allowed, and all values in 
#' this matrix should be positive or zero. For DNAm data, values should be 
#' beta-values.
#' 
#' @param method
#' Chioce of a reference-based method ('RPC','CBS','CP')
#' 
#' @param maxit
#' Only used in RPC mode, the limit of the number of IWLS iterations
#' 
#' @param nu.v
#' Only used in CBS mode. It is a vector of several candidate nu values. nu is 
#' parameter needed for nu-classification, nu-regression, and 
#' one-classification in svm. The best estimation results among all candidate nu 
#' will be automatically returned.
#' 
#' @param constraint
#' Only used in CP mode, you can choose either of 'inequality' or 'equality' 
#' normalization constraint. The default is 'inequality' (i.e sum of weights 
#' adds to a number less or equal than 1), which was implemented in 
#' Houseman et al (2012).
#' 
#' @return CP-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' ref: the reference centroid matrix used; dataREF: the subset of the input 
#' data matrix with only the probes defined in the reference matrix.
#' 
#' @return CBS-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' nu: a vector of 'best' nu-parameter for each sample; 
#' ref: the reference centroid matrix used;
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
#' 
#' @return RPC-mode
#' A list with the following entries: estF: a matrix of the estimated fractions;
#'  ref: the reference centroid matrix used; 
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
#' 
#' @references 
#' Teschendorff AE, Breeze CE, Zheng SC, Beck S. 
#' \emph{A comparison of reference-based algorithms for correcting cell-type 
#' heterogeneity in Epigenome-Wide Association Studies.}
#' BMC Bioinformatics (2017) 18: 105.
#' doi:\href{https://doi.org/10.1186/s12859-017-1511-5}{
#' 10.1186/s12859-017-1511-5}.
#' 
#' Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, 
#' Nelson HH, Wiencke JK, Kelsey KT. 
#' \emph{DNA methylation arrays as surrogate measures of cell mixture 
#' distribution.} 
#' BMC Bioinformatics (2012) 13: 86.
#' doi:\href{https://doi.org/10.1186/1471-2105-13-86}{10.1186/1471-2105-13-86}.
#' 
#' Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, Hoang CD, Diehn M, 
#' Alizadeh AA. 
#' \emph{Robust enumeration of cell subsets from tissue expression profiles.}
#' Nat Methods (2015) 12: 453-457.
#' doi:\href{https://doi.org/10.1038/nmeth.3337}{10.1038/nmeth.3337}.
#' 
#' @examples 
#' data(centDHSbloodDMC.m)
#' data(DummyBeta.m)
#' out.l <- epidish(DummyBeta.m, centDHSbloodDMC.m[,1:6], method = 'RPC')
#' frac.m <- out.l$estF
#' 
#' 
#' @export
#'     
epidish <- function(beta.m, ref.m, method = c("RPC", "CBS", "CP"), maxit = 50, nu.v = c(0.25, 
    0.5, 0.75), constraint = c("inequality", "equality")) {
    method <- match.arg(method)
    constraint <- match.arg(constraint)
    if (!method %in% c("RPC", "CBS", "CP")) 
        stop("Input a valid method!")
    if (method == "RPC") {
        out.o <- DoRPC(beta.m, ref.m, maxit)
    } else if (method == "CBS") {
        out.o <- DoCBS(beta.m, ref.m, nu.v)
    } else if (method == "CP") {
        if (!constraint %in% c("inequality", "equality")) {
            # make sure constraint must be inequality or equality
            stop("constraint must be inequality or equality when using CP!")
        } else out.o <- DoCP(beta.m, ref.m, constraint)
    }
    return(out.o)
}

### Reference-based methods

#' @importFrom MASS rlm
### RPC
DoRPC <- function(beta.m, ref.m, maxit) {
    map.idx <- match(rownames(ref.m), rownames(beta.m))
    rep.idx <- which(is.na(map.idx) == FALSE)
    data2.m <- beta.m[map.idx[rep.idx], , drop=FALSE]
    ref2.m <- ref.m[rep.idx, ]
    est.m <- matrix(nrow = ncol(data2.m), ncol = ncol(ref2.m))
    colnames(est.m) <- colnames(ref2.m)
    rownames(est.m) <- colnames(data2.m)
    for (s in seq_len(ncol(data2.m))) {
        rlm.o <- rlm(data2.m[, s] ~ ref2.m, maxit = maxit)
        coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m) + 1), 1]
        coef.v[which(coef.v < 0)] <- 0
        total <- sum(coef.v)
        coef.v <- coef.v/total
        est.m[s, ] <- coef.v
    }
    return(list(estF = est.m, ref = ref2.m, dataREF = data2.m))
}

#' @importFrom e1071 svm
### CIBERSORT
DoCBS <- function(beta.m, ref.m, nu.v) {
    map.idx <- match(rownames(ref.m), rownames(beta.m))
    rep.idx <- which(is.na(map.idx) == FALSE)
    
    data2.m <- beta.m[map.idx[rep.idx], , drop=FALSE]
    ref2.m <- ref.m[rep.idx, ]
    
    est.lm <- list()
    nui <- 1
    for (nu in nu.v) {
        est.m <- matrix(nrow = ncol(data2.m), ncol = ncol(ref2.m))
        colnames(est.m) <- colnames(ref2.m)
        rownames(est.m) <- colnames(data2.m)
        for (s in seq_len(ncol(data2.m))) {
            svm.o <- svm(x = ref2.m, y = data2.m[, s], scale = TRUE, type = "nu-regression", 
                kernel = "linear", nu = nu)
            coef.v <- t(svm.o$coefs) %*% svm.o$SV
            coef.v[which(coef.v < 0)] <- 0
            total <- sum(coef.v)
            coef.v <- coef.v/total
            est.m[s, ] <- coef.v
        }
        est.lm[[nui]] <- est.m
        nui <- nui + 1
    }
    
    #### select best nu
    rmse.m <- matrix(NA, nrow = ncol(beta.m), ncol = length(nu.v))
    for (nui in seq_along(nu.v)) {
        reconst.m <- ref2.m %*% t(est.lm[[nui]])
        s <- seq_len(ncol(beta.m))
        rmse.m[s, nui] <- sqrt(colMeans((data2.m[, s, drop = FALSE] - reconst.m[, s,  drop = FALSE])^2))
        message(nui)
    }
    colnames(rmse.m) <- nu.v
    nu.idx <- apply(rmse.m, 1, which.min)
    estF.m <- est.m
    for (s in seq_len(nrow(estF.m))) {
        estF.m[s, ] <- est.lm[[nu.idx[s]]][s, ]
    }
    return(list(estF = estF.m, nu = nu.v[nu.idx], ref = ref2.m, dataREF = data2.m))
}

#' @importFrom quadprog solve.QP   
### Houseman CP
DoCP <- function(beta.m, ref.m, constraint) {
    ### define D matrix
    nCT <- ncol(ref.m)
    D <- 2 * apply(ref.m, 2, function(x) colSums(x * ref.m))
    
    ### for inequality and equality, coe.v is different
    if (constraint == "inequality") {
        coe.v <- c(-1, 0)
    } else coe.v <- c(1, 1)
    
    ### define constraints
    A.m <- matrix(0, nrow = nCT, ncol = nCT)
    diag(A.m) <- rep(1, nCT)
    A.m <- cbind(rep(coe.v[1], nCT), A.m)
    b0.v <- c(coe.v[1], rep(0, nCT))
    
    ### define d-vector and solve for each sample
    nS <- ncol(beta.m)
    westQP.m <- matrix(NA, ncol = ncol(ref.m), nrow = nS)
    colnames(westQP.m) <- colnames(ref.m)
    rownames(westQP.m) <- colnames(beta.m)
    
    map.idx <- match(rownames(ref.m), rownames(beta.m))
    rep.idx <- which(is.na(map.idx) == FALSE)
    for (s in seq_len(nS)) {
        tmp.v <- beta.m[, s]
        d.v <- as.vector(2 * matrix(tmp.v[map.idx[rep.idx]], nrow = 1) %*% ref.m[rep.idx, 
            ])
        qp.o <- solve.QP(D, d.v, A.m, b0.v, meq = coe.v[2])
        westQP.m[s, ] <- qp.o$sol
        message(s)
    }
    
    return(list(estF = westQP.m, ref = ref.m[rep.idx, ], dataREF = beta.m[map.idx[rep.idx], 
        ]))
}
