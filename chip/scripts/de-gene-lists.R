source('scripts/de.R')

writeContrast <- function (type, contrast, name)
    write.table(contrast[order(contrast$padj), ],
                sprintf('results/de-%s/%s.tsv', type, sub('/', '-vs-', name)),
                quote = FALSE, sep = '\t', col.names = NA)

writeAllGeneLists <- function ()
    invisible(mapply(lp(writeContrast, 'genes'), trnaDeGenes, names(trnaDeGenes)))

if (! interactive()) {
    cat('# Generating gene list files\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    mkdir('results/de-genes')
    writeAllGeneLists()
}
