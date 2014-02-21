source('scripts/correlation.R')

writeAbundance <- function (data, name) {
    # Reorder columns
    cols <- paste(rep(tissues, each = length(stages)), stages, sep = '-')
    data <- data[, cols]
    write.table(data, file.path('results', name, ext = 'tsv'), sep = '\t',
                col.names = NA, quote = FALSE)
}

if (! interactive()) {
    cat('# Generating codon usage tables\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()
    trnaGroupFamilyAndType()

    writeAbundance(trnaByAcceptor, 'isoacceptors')
    writeAbundance(trnaByType[rownames(trnaByType) != 'SeCe', ], 'isotypes')
}
