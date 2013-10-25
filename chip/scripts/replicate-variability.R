source('../common/scripts/basic.R')
source('scripts/load-data.R')

plotPairwiseReplicates <- function (condition, data) {
    xy <- rownames(subset(trnaMapping, Condition == condition))
    x <- xy[1]
    y <- xy[2]
    lim <- max(data[, c(x, y)])
    plot(data[[x]], data[[y]],
         pch = 20, cex = 0.7, col = colors[1],
         xlab = readable(x), ylab = readable(y), log = 'xy',
         xlim = c(1, lim), ylim = c(1, lim),
         main = paste('Replicates for', readable(condition)))
    abline(0, 1, col = 'gray', untf = TRUE)

    rho <- cor(data[[x]], data[[y]], method = 'spearman')
    r2 <- cor(data[[x]], data[[y]], method = 'pearson')
    par(xlog = FALSE, ylog = FALSE, usr = c(0, 1, 0, 1))
    text(1, 0, bquote(atop(' ' ~ rho == .(sprintf('%.2f', rho)),
                           R^2 == .(sprintf('%.2f', r2)))),
         adj = c(1.1, -0.1))
}

plotReplicateVariability <- function () {
    plotOne <- function (condition) {
        on.exit(dev.off())
        pdf(file.path('plots', 'replicates', condition, ext = '.pdf'),
            width = 4, height = 4, family = plotFamily)
        plotPairwiseReplicates(condition, trnaNormData)
    }
    invisible(lapply(colnames(trnaNormDataCond), plotOne))
}

if (! interactive()) {
    cat('# Generating replicate variability plots\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    mkdir('plots/replicates')
    plotReplicateVariability()
}
