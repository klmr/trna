# Output aggregated normalised count data table with annotations
source('../common/scripts/basic.R')
source('scripts/load-data.R')

trnaWriteAnnotatedData <- function () {
    annotatedNormCond <- cbind(trnaNormDataCond,
                               trnaAnnotation[, c('Iso', 'Start', 'End', 'Strand')])

    write.table(annotatedNormCond,
                col.names = NA,
                file = 'results/aggregated-trnas.tsv',
                sep = '\t', quote = FALSE)
}

trnaWriteFlankedAnnotation <- function () {
    # ...
    stop('Not implemented.')
}

if (! interactive()) {
    cat('# Generating aggregated tRNA annotation\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    mkdir('results')
    trnaWriteAnnotatedData()
}
