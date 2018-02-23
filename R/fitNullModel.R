
setGeneric("fitNullModel2", function(x, ...) standardGeneric("fitNullModel2"))

setMethod("fitNullModel2",
          "data.frame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, group.var)
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                           group.idx=desmat$group.idx, ...)
          })

setMethod("fitNullModel2",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   cov.mat = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, group.var, sample.id)

              # subset or re-order cov.mat if necessary
              if (!is.null(cov.mat)) {
                  if (!is.list(cov.mat)) {
                      cov.mat <- list(A=cov.mat)
                  }
                  if (!is.null(sample.id)) {
                      cov.mat <- lapply(cov.mat, .orderSamples,
                                           orig.ids=x$sample.id, new.ids=sample.id)
                  }
              }
              
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=cov.mat,
                           group.idx=desmat$group.idx, ...)
          })

setMethod("fitNullModel2",
          "SeqVarData",
          function(x, ...) {
              fitNullModel2(sampleData(x), ...)
          })


invNormNullModel <- function(null.model, cov.mat = NULL, ...) {

    # subset or re-order cov.mat if necessary
    if (!is.null(cov.mat)) {
        if (!is.list(cov.mat)) {
            cov.mat <- list(A=cov.mat)
        }
        cov.mat <- lapply(cov.mat, function(x) {
            .orderSamples(x, orig.ids=rownames(x), new.ids=null.model$sample.id)
        })
    }

    updateNullModOutcome(null.model, covMatList=cov.mat, ...)
}


.orderSamples <- function(cov.mat, orig.ids, new.ids) {
    if (!is.null(rownames(cov.mat)) & !is.null(colnames(cov.mat))) {
        stopifnot(identical(rownames(cov.mat), colnames(cov.mat)))
        stopifnot(all(new.ids %in% rownames(cov.mat)))
    } else if (!is.null(rownames(cov.mat))) {
        stopifnot(all(new.ids %in% rownames(cov.mat)))
        colnames(cov.mat) <- rownames(cov.mat)
    } else if (!is.null(colnames(cov.mat))) {
        stopifnot(all(new.ids %in% colnames(cov.mat)))
        rownames(cov.mat) <- colnames(cov.mat)
    } else {
        warning("no dimnames given for cov.mat; assuming order of samples matches data frame")
        dimnames(cov.mat) <- list(orig.ids, orig.ids)
    }
    orig.ids <- as.character(orig.ids)
    new.ids <- as.character(new.ids)
    if (identical(rownames(cov.mat), new.ids)) {
        return(cov.mat)
    } else {
        return(cov.mat[new.ids, new.ids])
    }
}
