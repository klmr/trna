source('scripts/de.R')

trnaPlotCountMatrix <- function () {
    mkdir('plots/de')
    categories <- list(genes = trnaDeCounts,
                       acc = trnaAccDe$counts,
                       type = trnaTypeDe$counts)

    plotSingle <- function (data, name) {
        on.exit(dev.off())
        pdf(sprintf('plots/de/%s-counts.pdf', name), width = 7, height = 6,
            family = plotFamily)
        plotCountMatrix(data, 'Number of differentially expressed genes')
    }
    mapply(plotSingle, categories, names(categories))
}

if (! interactive()) {
    cat('# Generating heatmap of number of DE tRNA genes\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    trnaDetailedDe()
    trnaPlotCountMatrix()
}
