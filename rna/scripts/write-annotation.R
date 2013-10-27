source('../common/scripts/basic.R')
source('scripts/load-data.R')

mrnaWriteReplicates <- function () {
    annotatedNormData <- cbind(mrnaNormData,
                               mrnaAnnotation)
    write.table(annotatedNormData,
                col.names = NA,
                file = 'results/normalised-mrnas.tsv',
                sep = '\t', quote = FALSE)
}

mrnaWriteConditions <- function () {
    annotatedNormCond <- cbind(mrnaNormDataCond,
                               mrnaAnnotation)
    write.table(annotatedNormCond,
                col.names = NA,
                file = 'results/aggregated-mrnas.tsv',
                sep = '\t', quote = FALSE)
}

if (! interactive()) {
    cat('# Writing annotation for protein-coding genes\n')
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaNormalizeData()

    mkdir('results')
    mrnaWriteReplicates()
    mrnaWriteConditions()
}
