
setGeneric("fitNullModel2", function(x, ...) standardGeneric("fitNullModel2"))

setMethod("fitNullModel2",
          "AnnotatedDataFrame",
          function(x, outcome,
                   covars = NULL,
                   covMatList = NULL,
                   group.var = NULL,
                   sample.id = NULL,
                   ...) {
              desmat <- createDesignMatrix2(x, outcome, covars, group.var, sample.id)
              fitNullModel(y=desmat$y, X=desmat$X, covMatList=covMatList,
                           group.idx=desmat$group.idx, ...)
          })

setMethod("fitNullModel2",
          "SeqVarData",
          function(x, ...) {
              fitNullModel2(sampleData(x), ...)
          })
