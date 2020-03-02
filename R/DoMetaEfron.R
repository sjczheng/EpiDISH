#' @title An R-function to perform a meta-analysis over multiple studies using 
#' an empirical Bayes procedure by Efron followed by Stouffer method.
#' 
#' @aliases DoMetaEfron
#' 
#' @description
#' An R-function to perform a meta-analysis over multiple studies using 
#' an empirical Bayes procedure by Efron followed by Stouffer method.
#' 
#' @param stat.m
#' A matrix of signed statistics (e.g. t-statistics) with rows labeling genomic 
#' features (e.g. CpGs or genes) and columns labeling studies. rownames must be 
#' provided. 
#' 
#' @param pval.m
#' A matrix matched to stat.m containing the associated P-values, with rows 
#' labeling genomic features (e.g. CpGs or genes) and columns labeling studies.
#' 
#' @param bre
#' The number of breakpoints to divide statistics per study into bins. By 
#' default this is 120. See input argument for locfdr function from locfdr 
#' package.
#' 
#' @param df 
#' The number of degrees of freedom for fitting spline. By default this is 
#' 15. See input argument for locfdr function from locfdr package.
#' 
#' @param pct0
#' Percentage of statistics to use for fitting null. By default this is 0.25 
#' (i.e. 25\%).
#' 
#' @param plotlocfdr
#' Determines whether to plot output or not. By default this is set to 0 meaning
#'  no plot. See input argument for locfdr function from locfdr package.
#' 
#' @return meta.m
#' A matrix with rows as in stat.m, and with 3 columns labeling Stouffer's 
#' z-statistic, P-value and Benjamini-Hochberg adjusted P-value.
#' 
#' @import locfdr
#' 
#' @export
#' 

DoMetaEfron <- function(stat.m, pval.m, bre = 120, df = 15, pct0 = 0.25, plotlocfdr = 0){
  
  ### first transform p-values and signed t-stats to z-scores
  zMeta.m <- stat.m;
  zMeta.m[stat.m>0] <- qnorm(pval.m[stat.m>0]/2,lower.tail=FALSE);
  zMeta.m[stat.m<0] <- qnorm(pval.m[stat.m<0]/2,lower.tail=TRUE);
  
  ### for each study estimate z-stats of null
  zmodMeta.m <- zMeta.m;
  for(st in 1:ncol(zMeta.m)){
    locfdr.o <- locfdr(zMeta.m[,st], bre=bre, df=df, pct=0, pct0=pct0, nulltype=1, type=0, plot = plotlocfdr, sw=0);
    mu <- locfdr.o$fp0[3,1];
    sd <- locfdr.o$fp0[3,2];
    zmodMeta.m[,st] <- (zMeta.m[,st] - mu)/sd;
  }
  
  ### stouffer meta
  zMeta.v <- apply(zmodMeta.m,1,sum)/sqrt(ncol(stat.m));
  pMeta.v <- 2*(1-pnorm(abs(zMeta.v)));
  padjMeta.v <- p.adjust(pMeta.v,"BH");
  
  meta.m <- cbind(zMeta.v,pMeta.v,padjMeta.v);
  rownames(meta.m) <- rownames(stat.m);
  colnames(meta.m) <- c("z(Stouffer)","P(Stouffer)","AdjP(BH)");
  
  return(meta.m);
} ### EOF