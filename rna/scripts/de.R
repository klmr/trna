source('../common/scripts/basic.R')
source('scripts/load-data.R')

mrnaPairwiseDifferentialExpression <- function () {
    if (exists('mrnaDeGenes'))
        return()

    deGenesFile <- '../common/cache/de-genes.RData'
    mkdir('../common/cache')

    suppressWarnings(loaded <- tryCatch(load(deGenesFile), error = .('')))

    if (! identical(loaded, 'mrnaDeData')) {
        cat('Cache file not found -- re-generating DE gene lists.\n')

        mrnaDeData <- pairwiseDifferentialExpression(mrnaRawCounts, mrnaMapping, 0.05)
        save(mrnaDeData, file = deGenesFile)
    }

    mrnaDeResults <<- mrnaDeData$results
    mrnaDeGenes <<- mrnaDeData$de
    mrnaDeCount <<- mrnaDeData$counts
}
