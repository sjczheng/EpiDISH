---
title: "EpiDISH - Epigenetic Dissection of Intra-Sample-Heterogeneity - R package"
author: "Andrew E. Teschendorff, Shijie C. Zheng"
date: "`r Sys.Date()`"
package: "`r pkg_ver('EpiDISH')`"
output:
  BiocStyle::html_document
bibliography: EpiDISH.bib
vignette: >
  %\VignetteIndexEntry{Epigenetic Dissection of Intra-Sample-Heterogeneity - R package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

The **EpiDISH** package provides tools to infer the proportions of a priori known cell subtypes present in a sample representing a mixture of such cell-types. Inference proceeds via one of 3 methods (Robust Partial Correlations-RPC[@EpiDISH], Cibersort (CBS)[@CBS], Constrained Projection (CP)[@CP]), as determined by user.

For now, the package only includes one whole blood reference of 333 tsDHS-DMCs and 8 blood cell subtypes(B-cells, CD4+ T-cells, CD8+ T-cells, NK-cells, Monocytes, Neutrophils, Eosinophils, and Granulocytes. Note that Granulocytes consist of Neutrophils and Eosinophils.) described in [@EpiDISH]. This referecen dataset was based on 450k DNAm array; however, it could be directly used on both of 450k and EPIC array data. This package is under development and will offer reference-based inference for different tissue types. We will also include more algorithms in the future.


# How to use **EpiDISH** package

Using **EpiDISH** is quite simple. Here we use a small Illumina HumanMethylation450 BeadChip blood dataset(n=2) on GEO as an example.

You can download the dataset with *getGEO* function in **GEOquery** package and extract the whole beta value matrix.

```{r download, eval=F, echo=T, message=FALSE, warning=FALSE}
require(GEOquery)
require(Biobase)
GSE80559 <- getGEO("GSE80559")
beta.m <- exprs(GSE80559[[1]])
```

To reduce the package size and running time, we randomly selected 1000 probes from the beta value matrix(we let 330 of the probes be overlapped with the blood reference we provide.). The resulted *DummyBeta.m* is stored in the package. 

We load **EpiDISH** package, beta value matrix, and the whole blood reference dataset.
```{r load, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
library(EpiDISH)
data(centDHSbloodDMC.m)
data(DummyBeta.m)
```

Notice that *centDHSbloodDMC.m* has 8 columns. Granulocytes consist of Neutrophils and Eosinophils.
So, we only want to inlcude 7 columns(i.e B-cells, CD4+ T-cells, CD8+ T-cells, NK-cells, Monocytes, Neutrophils and Eosinophils) or 6 columns(i.e B-cells, CD4+ T-cells, CD8+ T-cells, NK-cells, Monocytes and Granulocytes). We go ahead and use *epidish* function with *RPC* mode to infer the proportions.
```{r infer, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
ref.m <- centDHSbloodDMC.m[,1:6]
out.l <- epidish(DummyBeta.m, ref.m, method = "RPC") 
```

Then, we check the output list. *estF* is the estimated cell fraction matrix. *ref* is the reference centroid matrix used; and *dataREF* is the input data matrix over the probes defined in the reference matrix.
```{r check, eval=TRUE, echo=T, message=FALSE, warning=FALSE}
out.l$estF
dim(out.l$ref)
dim(out.l$dataREF)
```

In this case, 330 out of 333 probes in the input reference matrix can be found in the inquiry matrix. So the *ref* is a $330*6$ matrix, while *dataREF* is a $330*2$ matrix. In QC step, we might remove bad probes; consequently, not all probes in the reference can be found in inquiry data. By checking *ref* and *dataREF*, we can extract the probes used to infer the proportions. If most of the probes in the reference cannot be found, the estimated proportions might be compromised.


# More info about different methods
We compared CP and RPC in [@EpiDISH]. And we also have a review article[@review] which summarized all methods tackling cell heterogeneity for DNAm data. Refers to references section for more details.

# Sessioninfo

```{r sessionInfo, echo=FALSE}
sessionInfo()
```

# References



