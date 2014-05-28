source('scripts/de.R')

if (! interactive()) {
    cat('Compute list of genes expressed in each stage\n')
    trnaLoadData()
    trnaPairwiseDiffentialExpression()

    activeGenes <- map(getExpressedtRNAs %.% countsForCondition, unique(trnaMapping$Condition))
    base <- 'results/active-genes'
    mkdir(base)
    map(.(name, data = writeLines(data, file.path(base, name, ext = 'txt'))),
        names(activeGenes), activeGenes) %|% invisible
}
