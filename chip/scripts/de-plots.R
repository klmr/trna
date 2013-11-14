source('scripts/de.R')

trnaPlotCountMatrix <- function () {
    mkdir('plots/de')
    on.exit(dev.off())
    pdf('plots/de/counts.pdf', width = 7, height = 6, family = plotFamily)
    plotCountMatrix(trnaDeCounts, 'Number of differentially expressed genes')
}

if (! interactive()) {
    cat('# Generating heatmap of number of DE tRNA genes\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    trnaPlotCountMatrix()
}
