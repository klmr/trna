source('../common/scripts/basic.R')
source('scripts/load-data.R')

trnaPcaCreatePlots <- function () {
    correlated <- -cor(trnaNormData, method = 'spearman')
    trnaPC <- prcomp(correlated)

    heatmapLab <- sapply(trnaMapping[colnames(trnaNormData), 'Condition'], readable)
    pdf(file.path(output, 'trna-heatmap.pdf'), width = 6, height = 6, family = plotFamily)
    heatmap(correlated, symm = TRUE, labRow = heatmapLab, labCol = heatmapLab,
            col = contrastColors)
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
