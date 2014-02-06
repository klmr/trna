source('../common/scripts/basic.R')
source('scripts/load-data.R')
local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/usage.R')
})

trnaGroupFamilyAndType <- function () {
    trnaByAcceptor <<- groupby(trnaNormDataCond, trnaAnnotation$Acceptor, sum)
    trnaByType <<- groupby(trnaNormDataCond, trnaAnnotation$Type, sum)
}

resampleAcceptorAbundance <- function () {
    samples <- 100
    cds <- trnaGetCountDataSet(trnaUnfilteredRawCounts)
    data <- trnaMergeReplicates(as.data.frame(counts(cds, normalized = TRUE)))
    data <- data[, grep('liver', colnames(data))]
    # Reorder columns by stages since they are unordered
    data <- data[, vapply(stages, grep, numeric(1), colnames(data))]

    sampleColumn <- function (col)
        map(.(. = shuffleRows(data[, col, drop = FALSE])), 1 : samples) %|%
        lp(do.call, cbind)
    acceptorSampleMatrix <- setNames(map(sampleColumn, colnames(data)),
                                     colnames(data))

    acceptorSampleMatrix <<- acceptorSampleMatrix
}
