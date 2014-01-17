source('scripts/de.R')

countsForCondition <- function (conditions)
    trnaRawCounts[, colnames(trnaRawCounts) %in%
                    rownames(trnaMapping)[trnaMapping$Condition %in% conditions]]

if (! interactive()) {
    trnaLoadData()
    trnaPairwiseDiffentialExpression()

    activeGenes <- map(getExpressedtRNAs %.% countsForCondition, unique(trnaMapping$Condition))
    base <- 'results/active-genes'
    mkdir(base)
    map(.(name, data = writeLines(data, file.path(base, name, ext = 'txt'))),
        names(activeGenes), activeGenes) %|% invisible
}
