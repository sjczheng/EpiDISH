#' @title 
#' Epigenetic Dissection of Intra-Sample-Heterogeneity
#' 
#' @aliases epidish
#'  
#' @description 
#' A reference-based function to infer the proportions of a priori known cell subtypes present in a sample representing a mixture of such cell-types. Inference proceeds via one of 3 methods (Robust Partial Correlations-RPC, Cibersort (CBS), Constrained Projection (CP)), as determined by user.
#' 
#' @param ref.m
#' A matrix of reference "centroids", i.e. representative molecular profiles, for a number of cell subtypes. rows label molecular features (e.g. CpGs,...) and columns label the cell-type. IDs need to be provided as rownames and colnames, respectively. No missing values are allowed, and all values in this matrix should be positive or zero. For DNAm data, values should be beta-values.
#' 
#' @param avdata.m
#' A data matrix with rows labeling the molecular features (should use same ID as in cent.m) and columns labeling samples (e.g. primary tumour specimens). No missing values are allowed and all values should be positive or zero. In the case of DNA methylation, these are beta-values.
#' 
#' @param method
#' Chioce of a reference-based method ("RPC","CBS","CP")
#' 
#' @param maxit
#' Used in RPC mode, the limit on the number of IWLS iterations
#' 
#' @param nu.v
#' This is only used for CBS mode. It is a vector of several nv values. nu is parameter needed for nu-classification, nu-regression, and one-classification in svm
#' 
#' @param constraint
#' For CP mode, you can choose either of "inequality" or "equality" normalization constraint. The default is "inequality" (i.e sum of weights adds to a number less or equal than 1), which was implemented in Houseman et al (2012).
#' 
#' @return CP-mode
#' A list with the following entries: estF: the estimated cell fraction matrix; ref: the reference centroid matrix used; dataREF: the input data matrix over the probes defined in the reference matrix
#' 
#' @return CBS-mode
#' A list with the following entries: estF: the estimated cell fraction matrix; nu: a vector of "best" nu-parameter for each sample; 
#' ref: the reference centroid matrix used;
#' dataREF: the input data matrix over the probes defined in the reference matrix
#' 
#' @return RPC-mode
#' ### A list with the following entries: estF: the estimated cell fraction matrix;
#' ref: the reference centroid matrix used; 
#' dataREF: the input data matrix over the probes defined in the reference matrix
#' 
#' @references 
#' Teschendorff AE, Breeze CE, Zheng SC, Beck S. 
#' \emph{A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies.}
#' BMC Bioinformatics (2017) 18: 105.
#' doi:\href{https://doi.org/10.1186/s12859-017-1511-5}{10.1186/s12859-017-1511-5}.
#' 
#' Houseman EA, Accomando WP, Koestler DC, Christensen BC, Marsit CJ, Nelson HH, Wiencke JK, Kelsey KT. 
#' \emph{DNA methylation arrays as surrogate measures of cell mixture distribution.} 
#' BMC Bioinformatics (2012) 13: 86.
#' doi:\href{https://doi.org/10.1186/1471-2105-13-86}{10.1186/1471-2105-13-86}.
#' 
#' Newman AM, Liu CL, Green MR, Gentles AJ, Feng W, Xu Y, Hoang CD, Diehn M, Alizadeh AA. 
#' \emph{Robust enumeration of cell subsets from tissue expression profiles.}
#' Nat Methods (2015) 12: 453-457.
#' doi:\href{https://doi.org/10.1038/nmeth.3337}{10.1038/nmeth.3337}.
#' 
#' @examples 
#'   #  library(EpiDISH)
#'   #  data(centDHSbloodDMC.m)
#'   #  out.l <- epidish(avdata.m, centDHSbloodDMC.m[,1:6], method = "RPC") 
#'   #  estF.m <- out.l$estF
#'   ## avdata.m is from samples you would like to infer weights of cell subtypes
#'   ## estF.m is the inferred fractions
#' @export
#' @import Biobase
#' 
epidish <- function(avdata.m,ref.m,method=c("RPC","CBS","CP"),maxit=50,nu.v=c(0.25,0.5,0.75),
                    constraint="inequality"){

    if(method=="RPC"){
       out.o <- DoRPC(avdata.m,ref.m,maxit);
    }
    else if (method=="CBS"){
       out.o <- DoCBS(avdata.m,ref.m,nu.v);
    }
    else if (method=="CP"){
       if (constraint %in% c("inequality", "equality")) {  # make sure constraint must be inequality or equality
         out.o <- DoCP(avdata.m,ref.m, constraint);
       } else stop("constraint must be inequality or equality!")
    }
    else {
        print("Input a valid method!");
    }

    return(out.o);
}
    
### Reference-based methods

#' @import MASS
### RPC
DoRPC <- function(avdata.m,ref.m,maxit){
    #require(MASS);
    map.idx <- match(rownames(ref.m),rownames(avdata.m));
    rep.idx <- which(is.na(map.idx)==FALSE);

    data2.m <- avdata.m[map.idx[rep.idx],];
    ref2.m <- ref.m[rep.idx,];

    est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
    colnames(est.m) <- colnames(ref2.m);
    rownames(est.m) <- colnames(data2.m);
    for(s in 1:ncol(data2.m)){
      rlm.o <- rlm( data2.m[,s] ~ ref2.m ,maxit=maxit)
      coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m)+1),1];
      coef.v[which(coef.v<0)] <- 0;
      total <- sum(coef.v);
      coef.v <- coef.v/total;
      est.m[s,] <- coef.v;
    }
    return(list(estF=est.m,ref=ref2.m,dataREF=data2.m));
}


#' @import e1071
###  CIBERSORT
DoCBS <- function(avdata.m,ref.m,nu.v=c(0.25,0.5,0.75)){

    #require(e1071);
    map.idx <- match(rownames(ref.m),rownames(avdata.m));
    rep.idx <- which(is.na(map.idx)==FALSE);

    data2.m <- avdata.m[map.idx[rep.idx],];
    ref2.m <- ref.m[rep.idx,];

    est.lm <- list();
    nui <- 1;
    for(nu in nu.v){
     est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
     colnames(est.m) <- colnames(ref2.m);
     rownames(est.m) <- colnames(data2.m);
     for(s in 1:ncol(data2.m)){
      svm.o <- svm(x=ref2.m,y=data2.m[,s],scale = TRUE, type="nu-regression", kernel ="linear", nu = nu);
      coef.v <- t(svm.o$coefs) %*% svm.o$SV;
      coef.v[which(coef.v<0)] <- 0;
      total <- sum(coef.v);
      coef.v <- coef.v/total;
      est.m[s,] <- coef.v;
     }
     est.lm[[nui]] <- est.m;
     print(nui);
     nui <- nui+1;
    }

   #### select best nu
   rmse.m <- matrix(NA,nrow=ncol(avdata.m),ncol=length(nu.v));
   for(nui in 1:length(nu.v)){
       reconst.m <- ref2.m %*% t(est.lm[[nui]]);
       for(s in 1:ncol(avdata.m)){
         rmse.m[s,nui] <- sqrt(mean((data2.m[,s] - reconst.m[,s])^2));
       }
       print(nui);
   }
   colnames(rmse.m) <- nu.v;
   nu.idx <- apply(rmse.m,1,which.min);
   estF.m <- est.m;    
   for(s in 1:nrow(estF.m)){
    estF.m[s,] <- est.lm[[nu.idx[s]]][s,];
   }
   return(list(estF=estF.m,nu=nu.v[nu.idx],ref=ref2.m,dataREF=data2.m));
}

#' @import quadprog
### Houseman CP
DoCP <- function(avdata.m,ref.m, constraint){
#require(quadprog);
### define D matrix
nCT <- ncol(ref.m);
D <- matrix(NA,nrow=nCT,ncol=nCT);
for(j in 1:nCT){
 for(k in 1:nCT){
   D[j,k] <- 2*sum(ref.m[,j]*ref.m[,k])
 }
}

if (constraint == "inequality") {
  ### define constraints
  A.m <- matrix(0,nrow=nCT+1,ncol=nCT);
  A.m[1,] <- -1;
  for(i in 1:nCT){
    A.m[1+i,i] <- 1;
  }
  A.m <- t(A.m);
  b0.v <- c(-1,rep(0,nCT));
  
  ### define d-vector and solve for each sample
  nS <- ncol(avdata.m);
  westQP.m <- matrix(NA,ncol=ncol(ref.m),nrow=nS);
  colnames(westQP.m) <- colnames(ref.m);
  rownames(westQP.m) <- colnames(avdata.m);
  
  match(rownames(ref.m),rownames(avdata.m)) -> map.idx;
  rep.idx <- which(is.na(map.idx)==FALSE);
  for(s in 1:nS){
    tmp.v <- avdata.m[,s];
    d.v <- as.vector(2*matrix(tmp.v[map.idx[rep.idx]],nrow=1) %*% ref.m[rep.idx,]);
    qp.o <- solve.QP(D,d.v,A.m,b0.v,meq=0);
    westQP.m[s,] <- qp.o$sol;
    print(s);
  }
} else # when constraint is equality
  {  
  ### define constraints
  A.m <- matrix(0,nrow=nCT+1,ncol=nCT);
  A.m[1,] <- 1;
  for(i in 1:nCT){
    A.m[1+i,i] <- 1;
  }
  A.m <- t(A.m);
  b0.v <- c(1,rep(0,nCT));
  
  ### define d-vector and solve for each sample
  nS <- ncol(avdata.m);
  westQP.m <- matrix(NA,ncol=ncol(ref.m),nrow=nS);
  colnames(westQP.m) <- colnames(ref.m);
  rownames(westQP.m) <- colnames(avdata.m);
  
  match(rownames(ref.m),rownames(avdata.m)) -> map.idx;
  rep.idx <- which(is.na(map.idx)==FALSE);
  for(s in 1:nS){
    tmp.v <- avdata.m[,s];
    d.v <- as.vector(2*matrix(tmp.v[map.idx[rep.idx]],nrow=1) %*% ref.m[rep.idx,]);
    qp.o <- solve.QP(D,d.v,A.m,b0.v,meq=1);
    westQP.m[s,] <- qp.o$sol;
    print(s);
  }
}


return(list(estF=westQP.m,ref=ref.m[rep.idx,],dataREF=avdata.m[map.idx[rep.idx],]));

}
