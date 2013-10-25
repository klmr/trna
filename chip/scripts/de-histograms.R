source('scripts/de.R')

#' @TODO These plots are partly redundant and a bit chaotic. Find a better representation for this.

plotHistDe <- function (dataDE, path, name, ...) {
    on.exit(dev.off())
    pdf(sprintf('%s.pdf', path), family = plotFamily, ...)
    plotOne <- function (axis, axisName) {
        if (is.null(axis))
            return()
        hist(axis$pval, breaks = 100, col = last(colors),
             main = sprintf('Differential expression %s%s', name, axisName),
             ylim = c(0, 50), xlab = 'p-value')
        text(0.2, 10, sprintf('Significant: %d tRNAs', count(axis$padj < 0.05)))
    }
    invisible(mapply(plotOne, dataDE, names(dataDE)))
}

plotAllHistograms <- function () {
    filtered <- function (tissue, a, b) {
        if (missing(b))
            trnaDeResults[[tissue]][[a]]
        else # assume missing(a)
            cdict(lapply(trnaDeResults[[tissue]], item(b)),
                  trnaDeResults[[tissue]][[b]])
    }

    allStages <- list(liver = filtered('liver', 'e15.5'),
                      brain = filtered('brain', 'e15.5'))
    adultStages <- list(liver = filtered('liver', b = 'P0.5'),
                        brain = filtered('brain', b = 'P0.5'))

    plotHistDe(allStages$brain, 'plots/de-hist/brain_e15.5', 'brain ', width = 8, height = 8)
    plotHistDe(allStages$liver, 'plots/de-hist/liver_e15.5', 'liver ', width = 8, height = 8)
    plotHistDe(adultStages$brain, 'plots/de-hist/brain_P0.5', 'brain ', width = 8, height = 8)
    plotHistDe(adultStages$liver, 'plots/de-hist/liver_P0.5', 'liver ', width = 8, height = 8)
    plotHistDe(trnaTissueDeResults, 'plots/de-hist/stages_liver_brain', '', width = 8, height = 8)
}

if (! interactive()) {
    cat('# Generating DE histograms\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    trnaTissueDifferentialExpression()
    mkdir('plots/de-hist')
    plotAllHistograms()
}
