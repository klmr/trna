source('../common/scripts/basic.R')
source('scripts/load-data.R')

library(gplots)

trnaPcaCreatePlots <- function () {
    correlated <- cor(trnaNormData, method = 'spearman')
    trnaPC <- prcomp(t(correlated), scale. = TRUE)
    # Flip data points on their axes to make the plot have the same orientation
    # as the mRNA plot. This doesn't affect the interpretation of the PCA.
    trnaPC$x <- -trnaPC$x
    #trnaPC <- prcomp(t(trnaNormData), scale. = TRUE)

    heatmapLab <- sapply(trnaMapping[colnames(trnaNormData), 'Condition'], readable)
    pdf(file.path(output, 'trna-heatmap.pdf'), width = 6, height = 6, family = plotFamily)
    heatmap.2(correlated, symm = TRUE, Colv = TRUE, labRow = heatmapLab,
              labCol = heatmapLab, col = contrastColors, trace = 'none',
              density.info = 'none')
    dev.off()

    pdf(file.path(output, 'trna-pca.pdf'), width = 6, height = 6, family = plotFamily)
    plotProgression(trnaPC, 'Tissue', 'Stage', tissues, stages,
                    trnaMapping, '', tissueColor, legend = FALSE)
    dev.off()
}

if (! interactive()) {
    cat('# Generating correlation heatmap and PCA\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    output <- 'plots/distribution'
    mkdir(output)
    trnaPcaCreatePlots()
}
