source('../common/scripts/basic.R')
source('scripts/load-data.R')

createPlots <- function () {
    correlated <- cor(mrnaNormData, method = 'spearman')
    pc <- prcomp(correlated)

    #heatmapLab <- sprintf('%s %s', mrnaMapping[colnames(mrnaNormData), 'Tissue'], mrnaMapping[, 'Stage'])
    heatmapLab <- mrnaMapping[colnames(mrnaNormData), 'Condition']
    pdf(file.path(output, 'mrna-heatmap.pdf'), width = 6, height = 6)
    heatmap(correlated, symm = TRUE, labRow = heatmapLab, labCol = heatmapLab)
    dev.off()

    pdf(file.path(output, 'mrna-pca.pdf'), width = 6, height = 6)
    plotProgression(pc, 'Tissue', 'Stage', tissues, stages,
                    mrnaMapping, '', tissueColor, legend = FALSE)
    dev.off()
}

output <- 'plots/distribution'

if (! interactive()) {
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaNormalizeData()

    mkdir(output)
    createPlots()
}
