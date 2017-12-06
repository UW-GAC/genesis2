setGeneric("assocTestSeq2", function(gdsobj, ...) standardGeneric("assocTestSeq2"))

## should allele frequency filter be an argument, or set ahead of time?
## if an argument, can make sure freq is calcuated on sample set used
## compare seqSetFilterCond to calcluating MAF separately
## previously, returned rows in variantInfo for high-freq variants but set weight to 0
##  - why do this instead of just filtering?

setMethod("assocTestSeq2",
          "SeqVarIterator",
          function(gdsobj, nullModel, AF.max=1,
                   weight.beta=c(1,1), weight.user=NULL,
                   test=c("Burden", "SKAT"),
                   burden.test=c("Score", "Wald"), rho=0,
                   pval.method=c("davies", "kuonen", "liu"),
                   verbose=TRUE) {

              # check argument values
              test <- match.arg(test)
              burden.test <- match.arg(burden.test)
              pval.method <- match.arg(pval.method)

              # filter samples to match null model
              .setFilterNullModel(gdsobj, nullModel, verbose=verbose)
              
              nullPrep <- nullModelTestPrep(nullModel)
              
              # results
              res <- list()
              res.var <- list()
              i <- 1
              iterate <- TRUE
              while (iterate) {
                  var.info <- variantInfo(gdsobj, alleles=FALSE, expanded=TRUE)
                  
                  geno <- expandedAltDosage(gdsobj, use.names=FALSE)
                  
                  # allele frequency
                  freq <- .alleleFreq(gdsobj, geno)

                  # exclude monomorphic variants
                  mono <- freq %in% c(0,1)
                  # exclude variants with freq > max
                  excl <-  mono | freq > AF.max
                  if (any(excl)) {
                      var.info <- var.info[!excl,,drop=FALSE]
                      geno <- geno[,!excl,drop=FALSE]
                      freq <- freq[!excl]
                  }

                  # number of variant sites
                  n.site <- length(unique(var.info$variant.id))

                  # number of alternate alleles
                  n.alt <- sum(geno, na.rm=TRUE)
                  
                  # number of samples with observed alternate alleles > 0
                  n.sample.alt <- sum(rowSums(geno, na.rm=TRUE) > 0)
               
                  # number of non-missing samples
                  n.obs <- colSums(!is.na(geno))
                  
                  # weights
                  if (is.null(weight.user)) {
                      # Beta weights
                      weight <- .weightFromFreq(freq, weight.beta)
                  } else {
                      # user supplied weights
                      weight <- variantData(gdsobj)[[weight.user]][expandedVariantIndex(gdsobj)]
                      weight <- weight[!excl]
                  }
                  
                  res[[i]] <- data.frame(n.site, n.alt, n.sample.alt)
                  res.var[[i]] <- cbind(var.info, n.obs, freq, weight)
                  
                  if (n.site > 0) {
                      # mean impute missing values
                      if (any(n.obs < nrow(geno))) {
                          geno <- .meanImpute(geno, freq)
                      }

                      # do the test
                      assoc <- testVariantSet(nullPrep, G=geno, weights=weight, test=test, 
                                              burden.test=burden.test, rho=rho,
                                              pval.method=pval.method)
                      res[[i]] <- cbind(res[[i]], assoc)
                  }

                  i <- i + 1
                  iterate <- iterateFilter(gdsobj, verbose=verbose)
              }

              list(results=bind_rows(res), variantInfo=res.var)
          })
