
setGeneric("fitNullModel2", function(x, ...) standardGeneric("fitNullModel2"))

setMethod("fitNullModel2",
          "data.frame",
          function(x, outcome,
                   covars = NULL,
                   covMatList = NULL,
                   group.var = NULL,
                   ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, group.var)
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=covMatList,
                           group.idx=desmat$group.idx, ...)
          })

setMethod("fitNullModel2",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   covMatList = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, group.var, sample.id)

              # subset or re-order covMatList if necessary
              if (!is.null(covMatList)) {
                  if (!is.list(covMatList)) {
                      covMatList <- list(A=covMatList)
                  }
                  if (!is.null(sample.id)) {
                      covMatList <- lapply(covMatList, .orderSamples,
                                           orig.ids=x$sample.id, new.ids=sample.id)
                  }
              }
              
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=covMatList,
                           group.idx=desmat$group.idx, ...)
          })

setMethod("fitNullModel2",
          "SeqVarData",
          function(x, ...) {
              fitNullModel2(sampleData(x), ...)
          })

.orderSamples <- function(covMatList, orig.ids, new.ids) {
    if (!is.null(rownames(covMatList)) & !is.null(colnames(covMatList))) {
        stopifnot(identical(rownames(covMatList), colnames(covMatList)))
        stopifnot(all(new.ids %in% rownames(covMatList)))
    } else if (!is.null(rownames(covMatList))) {
        stopifnot(all(new.ids %in% rownames(covMatList)))
        colnames(covMatList) <- rownames(covMatList)
    } else if (!is.null(colnames(covMatList))) {
        stopifnot(all(new.ids %in% colnames(covMatList)))
        rownames(covMatList) <- colnames(covMatList)
    } else {
        warning("no dimnames given for covMatList; assuming order of samples matches data frame")
        dimnames(covMatList) <- list(orig.ids, orig.ids)
    }
    orig.ids <- as.character(orig.ids)
    new.ids <- as.character(new.ids)
    if (identical(rownames(covMatList), new.ids)) {
        return(covMatList)
    } else {
        return(covMatList[new.ids, new.ids])
    }
}
