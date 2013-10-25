source('scripts/de.R')

writeAllGeneLists <- function () {
    writeContrast <- function (contrast, name)
        write.table(contrast[order(contrast$padj), ],
                    sprintf('results/de-genes/%s.tsv', sub('/', '-vs-', name)),
                    quote = FALSE, sep = '\t', col.names = NA)

    invisible(mapply(writeContrast, trnaDeGenes, names(trnaDeGenes)))
}

if (! interactive()) {
    cat('# Generating gene list files\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    mkdir('results/de-genes')
    writeAllGeneLists()
}
