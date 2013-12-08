source('scripts/usage.R')

plotRnaCorrelation <- function () {
    methods <- list(gene = mrnaNormDataCond,
                    codon = codonUsageData,
                    aa = aaUsageData)
    titles <- list(gene = 'Protein-coding gene expression correlation',
                   codon = 'Codon usage correlation',
                   aa = 'Amino acid usage correlation')
    generateCorrelationPlot <- function (name, data, title) {
        on.exit(dev.off())
        pdf(file.path('plots', 'correlation',
                      sprintf('%s-usage', name), ext = 'pdf'))
        plotCorrelationsFor(data, title)
    }

    invisible(mapply(generateCorrelationPlot, names(methods), methods, titles))
}

writeCodonUsageData <- function()
    write.table(codonUsageData, file = 'results/codon-usage.dat')

plotSpiderWeb <- function (type, codonUsageData, aaUsageData) {
    relativeData <- function (data) {
        sums <- apply(data, COLS, sum)
        t(apply(data, ROWS, function (row) row / sums))
    }

    plotRadial <- function (data, labels, main, ...) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 0.2))
        radial.plot(t(data), labels = labels,
                    rp.type = '', show.grid.labels = 3, main = main, ...)
        lwd <- c(rep(2, ncol(data) - 1), 5)
        radial.plot(t(data),
                    rp.type = 'p', line.col = tcolors, lwd = lwd,
                    show.grid = FALSE, show.radial.grid = FALSE, add = TRUE, ...)
        plot.new()
        legend('center', legend = c(stages, 'back­\nground'), fill = tcolors,
               border = NA, bty = 'n', xpd = NA)
    }

    require(plotrix)
    tcolors <- c(colors[indices(stages)], '#00000080')

    for (tissue in tissues) {
        data <- aaUsageData[, grep(tissue, colnames(aaUsageData))]
        data$Background <- aaBackgroundDist$Count
        # Enforce uniform oder between tRNA and mRNA plots.
        data <- data[aminoAcids$Short, ]
        rownames(data) <- aminoAcids[aminoAcids$Short == rownames(data), 'Long']
        relative <- relativeData(data)
        pdf(sprintf('plots/usage/%s%s.pdf', type, tissue),
            width = 6, height = 6, family = plotFamily)
        plotRadial(relative, rownames(data),
                   radial.lim = c(0, 0.1),
                   main = sprintf('Amino acid usage across stages in %s\n', tissue))
        dev.off()

        # Codons of Arginine
        rcodons <- rownames(subset(geneticCode, AA == 'R'))
        data <- codonUsageData[rcodons, grep(tissue, colnames(codonUsageData))]
        data$Background <- codonBackgroundDist[rcodons, 'Count']
        relative <- relativeData(data)
        pdf(sprintf('plots/usage/%sarginine-%s.pdf', type, tissue),
            width = 6, height = 6, family = plotFamily)
        plotRadial(relative, rownames(data), radial.lim = c(0, 0.3),
                   main = sprintf('Arginine codon usage across stages in %s\n', tissue))
        dev.off()
    }
}

if (! interactive()) {
    cat('# Generating mRNA codon usage data\n')
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaNormalizeData()
    generateCodonUsageData()
    loadAminoAcids()
    loadCodonMap()
    generateCodonBackgroundDist()

    writeCodonUsageData()

    cat('# Generating mRNA expression & usage correlation plots\n')
    mkdir('plots/correlation')
    plotRnaCorrelation()

    cat('# Generating mRNA codon usage plots\n')
    mkdir('plots/usage')
    plotSpiderWeb('', codonUsageData, aaUsageData)
}
