## ----load, eval=TRUE, echo=T, message=FALSE, warning=FALSE-----------------
library(EpiDISH)
data(centEpiFibIC.m)
data(DummyBeta.m)

## ----infer, eval=TRUE, echo=T, message=FALSE, warning=FALSE----------------
out.l <- epidish(beta.m = DummyBeta.m, ref.m = centEpiFibIC.m, method = "RPC") 

## ----check, eval=TRUE, echo=T, message=FALSE, warning=FALSE----------------
out.l$estF
dim(out.l$ref)
dim(out.l$dataREF)

## ----hepidish, eval=TRUE, echo=T, message=FALSE, warning=FALSE-------------
data(centBloodSub.m)
frac.m <- hepidish(beta.m = DummyBeta.m, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m[,c(1, 2, 5)], h.CT.idx = 3, method = 'RPC')
frac.m

## ----celldmc, eval=TRUE, echo=T, message=FALSE, warning=FALSE--------------
pheno.v <- rep(c(0, 1), each = 5)
celldmc.o <- CellDMC(DummyBeta.m, pheno.v, frac.m)

## ----dmct, eval=TRUE, echo=T, message=FALSE, warning=FALSE-----------------
head(celldmc.o$dmct)

## ----sessionInfo, echo=FALSE-----------------------------------------------
sessionInfo()

