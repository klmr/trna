source('scripts/de.R')

mrnaWriteDeGenes <- function () {
    for (tissue in tissues) {
        base <- file.path('results/de', tissue)
        for (a in names(mrnaDeGenes[[tissue]])) {
            contrastA <- mrnaDeGenes[[tissue]][[a]]
            for (b in names(contrastA)) {
                if (is.null(contrastA[[b]])) {
                    cat('Nothing for', tissue, a, b, '\n')
                    next
                }
                write.table(contrastA[[b]], sep = '\t', quote = FALSE,
                            file = file.path(base, sprintf('%s-%s', a, b),
                                             ext = 'tsv'))
            }
        }
    }
}

if (! interactive()) {
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaPairwiseDifferentialExpression()
    mrnaWriteDeGenes()
}
