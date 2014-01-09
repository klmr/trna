source('../common/scripts/basic.R')
source('scripts/load-data.R')

mrnaPcaCreatePlots <- function () {
    correlated <- cor(mrnaNormData, method = 'spearman')
    pc <- prcomp(correlated)

    heatmapLab <- sapply(mrnaMapping[colnames(mrnaNormData), 'Condition'], readable)
    pdf(file.path(output, 'mrna-heatmap.pdf'), width = 6, height = 6)
    heatmap(correlated, symm = TRUE, labRow = heatmapLab, labCol = heatmapLab,
            col = contrastColors)
    dev.off()

    pdf(file.path(output, 'mrna-pca.pdf'), width = 6, height = 6)
    plotProgression(pc, 'Tissue', 'Stage', tissues, stages,
                    mrnaMapping, '', tissueColor, legend = FALSE)
    dev.off()
}

output <- 'plots/distribution'

if (! interactive()) {
    cat('# Generating heatmap and PCA of mRNA expression correlation\n')
    mrnaLoadData()
    mrnaNormalizeData()

    mkdir(output)
    mrnaPcaCreatePlots()
}
