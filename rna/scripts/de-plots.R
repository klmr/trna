source('scripts/de.R')

mrnaPlotCountMatrix <- function () {
    mkdir('plots/de')
    on.exit(dev.off())
    pdf('plots/de/counts.pdf', width = 7, height = 6, family = plotFamily)
    plotCountMatrix(mrnaDeCount, 'Number of differentially expressed genes')
}

if (! interactive()) {
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaPairwiseDifferentialExpression()
    mrnaPlotCountMatrix()
}
