context("epidish")

test_that("# of columns in output should equal to that of reference", {
    
    data(centDHSbloodDMC.m)
    data(DummyBeta.m)
    idx <- sample(1:8, sample(2:8, 1), replace = F)
    ref.m <- centDHSbloodDMC.m[, idx]
    estF.m <- epidish(DummyBeta.m, ref.m, method = "RPC")$estF
    
    dim1 <- ncol(ref.m)
    dim2 <- ncol(estF.m)
    expect_identical(dim1, dim2)
})
