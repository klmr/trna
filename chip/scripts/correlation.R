source('../common/scripts/basic.R')
source('scripts/load-data.R')
local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/correlation.R')
})

prepareForPairsPlot <- function (df) {
    mat <- list(brain = df[, grep('brain', colnames(df))],
                liver = df[, grep('liver', colnames(df))])
    newcols <- gsub('brain_', '', colnames(mat[[1]]))
    colnames(mat[[1]]) <- newcols
    colnames(mat[[2]]) <- newcols
    data.frame(do.call(rbind, mat[tissues]),
               type = rep(tissues, each = nrow(df)))
}

trnaGroupFamilyAndType <- function () {
    conditions <- unique(trnaMapping$Condition)
    trnaStageCount <<- prepareForPairsPlot(trnaNormDataCond)

    trnaByAcceptor <- aggregate(trnaNormDataCond,
                                by = list(Condition = trnaAnnotation[rownames(trnaNormDataCond), 'Acceptor']),
                                mean)
    rownames(trnaByAcceptor) <- trnaByAcceptor$Condition
    trnaByAcceptor$Condition <- NULL
    trnaByAcceptor <<- prepareForPairsPlot(trnaByAcceptor)

    trnaByType <- aggregate(trnaNormDataCond,
                            by = list(Condition = trnaAnnotation[rownames(trnaNormDataCond), 'Type']),
                            mean)
    rownames(trnaByType) <- trnaByType$Condition
    trnaByType$Condition <- NULL
    trnaByType <<- prepareForPairsPlot(trnaByType)
}

relativeData <- function (data) {
    sums <- apply(data, COLS, sum)
    t(apply(data, ROWS, function (row) row / sums))
}
