source('scripts/usage.R')

plotRnaCorrelation <- function () {
    methods <- list(gene = mrnaNormDataCond,
                    codon = codonUsageData,
                    aa = aaUsageData)
    titles <- list(gene = 'Protein­coding gene expression correlation',
                   codon = 'Codon usage correlation',
                   aa = 'Amino acid usage correlation')
    generateCorrelationPlot <- function (name, data, title) {
        on.exit(dev.off())
        pdf(file.path('plots', 'correlation',
                      sprintf('%s-usage', name), ext = 'pdf'))
        data <- map(.(tissue = cor(data[, grep(tissue, colnames(data))],
                                   method = 'spearman')), tissues)
        data <- map(p(`diag<-`, NA), data)
        data <- map(p(`colnames<-`, stages), data)
        data <- map(p(`rownames<-`, stages), data)
        plotCountMatrix(data, title, '%0.2f')
    }

    invisible(mapply(generateCorrelationPlot, names(methods), methods, titles))
}

writeCodonUsageData <- function()
    write.table(codonUsageData, file = 'results/codon-usage.dat')

plotSpiderWeb <- function (type, codonUsageData, aaUsageData) {
    plotRadial <- function (data, labels, main, ...) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 0.2))
        lwd <- c(rep(2, ncol(data) - 1), 5)
        radial.plot(data, labels = labels, main = main, line.col = tcolors,
                    lwd = lwd, show.grid.labels = 3, ...)
        plot.new()
        legend('center', legend = c(stages, 'back­\nground'), fill = tcolors,
               border = NA, bty = 'n', xpd = NA)
    }

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

plotCodonBackground <- function () {
    tcolors <- c(gray.colors(5), colors[1])
    lwd <- c(rep(2, 5), 3)
    local({
        on.exit(dev.off())
        pdf('plots/usage/aa-background.pdf', width = 6, height = 6,
            family = plotFamily)
        data <- relativeData(overallAaBackground)[aminoAcids$Short, ]
        data <- cbind(data[, -1], data[, 1])
        radial.plot(data, labels = rownames(data), line.col = tcolors,
                    show.grid.labels = 3, lwd = lwd,
                    main = 'Amino acid usage background for all reading frames')
    })

    local({
        on.exit(dev.off())
        pdf('plots/usage/codon-background.pdf', width = 10, height = 4,
            family = plotFamily)
        data <- as.data.frame(relativeData(overallCodonBackground))
        sortedCodons <- rownames(data)[order(geneticCode[rownames(data), 1])]
        groups <- unlist(map(length, split(sortedCodons, geneticCode[sortedCodons, 1])))
        xcolors <- unlist(map(rep, sapply(c(0.5, 0.8), gray), groups))
        ylim <- c(0, max(data))
        barpos <- barplot(data[sortedCodons, 1], col = xcolors, border = NA, xaxt = 'n',
                          ylim = ylim, ylab = 'Frequency',
                          main = 'Codon usage background for all reading frames')
        axispos <- groupby(barpos, geneticCode[sortedCodons, 1], mean)[[1]]
        axis(1, at = axispos, names(groups), tick = FALSE)
        map(.(x, col = lines(barpos[, 1], x, col = col)),
            data[sortedCodons, -1], colors[indices(data[, -1])])
    })
}

plotCodonUsageForAA <- function (aa, data, name) {
    codons <- rownames(subset(geneticCode, AA == aa))
    if (length(codons) == 1)
        return()

    long <- subset(aminoAcids, Short == aa)$Long
    main <- sprintf('Simulated %s tripled codon usage in %s',
                    long, readable(name))
    plotCodonBarplot(relativeData(data[codons, ]), main)
}

plotCodonSampling <- function () {
    map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/codons-%s.pdf', name))
        map(p(plotCodonUsageForAA, data, name), aminoAcids$Short)
    }), codonSampleMatrix, names(codonSampleMatrix)) %|% invisible
}

plotAAUsage <- function (data, name) {
    prepare <- function (data) {
        data <- groupby(data, geneticCode[rownames(data), 'AA'])
        # Enforce uniform oder between tRNA and mRNA plots.
        data <- data[aminoAcids$Short, ]
        rownames(data) <- aminoAcids[aminoAcids$Short == rownames(data), 'Long']
        relativeData(data)
    }

    n <- ncol(data[[1]])
    data <- do.call(cbind, map(prepare, data))

    radial.plot(data, labels = rownames(data),
                line.col = rep(transparent(colors, 0.2), each = n),
                lwd = 2, show.grid.labels = 3, radial.lim = c(0, 0.1),
                main = 'Amino acid usage with resampled expression')
}

plotAminoAcidSampling <- function () {
    on.exit(dev.off())
    pdf('plots/usage-sampling/amino-acids.pdf')
    plotAAUsage(codonSampleMatrix, readable(names(codonSampleMatrix)))
}

if (! interactive()) {
    cat('# Generating mRNA codon usage data\n')
    mrnaLoadData()
    mrnaNormalizeData()
    generateCodonUsageData()
    loadAminoAcids()
    loadCodonMap()
    generateStableCodonUsageData()
    generateCodonBackgroundDist()
    generateCodonBackgroundUsage()

    writeCodonUsageData()

    cat('# Generating mRNA expression & usage correlation plots\n')
    mkdir('plots/correlation')
    plotRnaCorrelation()

    cat('# Generating mRNA codon usage plots\n')
    mkdir('plots/usage')
    plotSpiderWeb('', codonUsageData, aaUsageData)

    cat('# Generate stable mRNA codon usage plots\n')
    plotSpiderWeb('stable-', stableCodonUsageData, stableAaUsageData)

    cat('# Generate background mRNA codon usage plots\n')
    plotCodonBackground()

    cat('# Generate shuffled expression codon profiles\n')
    resampleCodonUsage()
    mkdir('plots/usage-sampling')
    plotCodonSampling()
    plotAminoAcidSampling()
}
