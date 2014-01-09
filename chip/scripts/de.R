source('../common/scripts/basic.R')
source('../common/scripts/de.R')
source('scripts/load-data.R')

trnaPairwiseDiffentialExpression <- function () {
    if (exists('trnaDeGenes'))
        return()

    results <- pairwiseDifferentialExpression(trnaRawCounts, trnaMapping, 0.05)

    trnaDeResults <<- results$results
    trnaDeGenes <<- results$de
    trnaDeCounts <<- results$counts
}

trnaDetailedDe <- function () {
    if (exists('trnaAccDe'))
        return()

    getData <- function (x) groupby(trnaRawCounts, trnaAnnotation[[x]])

    data <- map(getData, c('Acceptor', 'Type'))
    trnaAccDe <<- pairwiseDifferentialExpression(data$Acceptor, trnaMapping, 0.05)
    trnaTypeDe <<- pairwiseDifferentialExpression(data$Type, trnaMapping, 0.05)
}

trnaTissueDifferentialExpression <- function () {
    if (exists('trnaTissueDeGenes'))
        return()

    doDe <- function (stage) {
        cond <- p(paste, stage, sep = '-')
        nbinomTest(trnaCds, cond('liver'), cond('brain'))
    }
    trnaTissueDeResults <<- lapply(stages, p(sorted, 'padj') %.% doDe)
}
