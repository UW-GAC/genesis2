context("single variant tests")
library(SeqVarTools)

test_that("assocTestMM2", {
    svd <- .testData()
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    assoc <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    seqResetFilter(svd, verbose=FALSE)
    expect_equal(unique(assoc$variant.id), seqGetData(svd, "variant.id"))
    seqClose(svd)
})

test_that("assocTestMM2 - sample selection", {
    svd <- .testData()
    samp <- sampleData(svd)$sample.id[sample(1:nrow(sampleData(svd)), 50)]
    iterator <- SeqVarBlockIterator(svd, variantBlock=500, verbose=FALSE)
    nullmod <- fitNullModel2(iterator, outcome="outcome", covars=c("sex", "age"), sample.id=samp, verbose=FALSE)
    expect_equal(nrow(nullmod$model.matrix), 50)
    assoc <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    expect_equal(max(assoc$n.obs), 50)
    seqClose(svd)
})


test_that("assocTestMM2 matches regression", {
    svd <- .testData()
    
    # multiallelic variants are handled differently
    snv <- isSNV(svd, biallelic=TRUE)
    seqSetFilter(svd, variant.sel=snv, verbose=FALSE)
    assoc1 <- regression(svd, outcome="outcome", covar=c("sex", "age"))
    
    nullmod <- fitNullModel2(svd, outcome="outcome", covars=c("sex", "age"), verbose=FALSE)
    iterator <- SeqVarBlockIterator(svd, verbose=FALSE)
    assoc2 <- assocTestMM2(iterator, nullmod, verbose=FALSE)
    expect_equal(nrow(assoc1), nrow(assoc2))
    expect_equal(assoc1$variant.id, assoc2$variant.id)
    expect_equal(assoc1$n, assoc2$n.obs)
    expect_equal(assoc1$freq, 1-assoc2$freq)
    ## this won't match exactly, because missing data is handled differently
    ## assocTestMM2 imputes to the mean, while regression drops missing data
    #expect_equal(assoc1$Est, -assoc2$Est)
    #expect_equal(assoc1$SE, assoc2$Est.SE)
    #expect_equal(assoc1$Wald.Stat, (assoc2$Wald.Stat)^2)
    #expect_equal(assoc1$Wald.Pval, assoc2$Wald.pval, tolerance=.1)
    
    seqClose(svd)
})
