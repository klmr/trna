source('scripts/de.R')

writeContrast <- function (type, contrast, name)
    write.table(contrast[order(contrast$padj), ],
                sprintf('results/de-%s/%s.tsv', type, sub('/', '-vs-', name)),
                quote = FALSE, sep = '\t', col.names = NA)

writeAllGeneLists <- function () {
    categories <- list(genes = trnaDeGenes,
                       acc = trnaAccDe$de,
                       type = trnaTypeDe$de)

    writeCategory <- function (category, name) {
        mkdir(sprintf('results/de-%s', name))
        mapply(lp(writeContrast, name), category, names(category))
    }

    invisible(mapply(writeCategory, categories, names(categories)))
}

if (! interactive()) {
    cat('# Generating gene list files\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    trnaDetailedDe()
    writeAllGeneLists()
}
