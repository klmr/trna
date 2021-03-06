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
        plotCountMatrix(data, title, '%0.2f', commonScale = globalScale)
    }

    #' @TODO Do not hard-code the value
    globalScale <- c(0.9, 1, 0.1)
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
        data$Background <- aaBackgroundDist[, 1]
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
        data$Background <- codonBackgroundDist[rcodons, 1]
        relative <- relativeData(data)
        pdf(sprintf('plots/usage/%sarginine-%s.pdf', type, tissue),
            width = 6, height = 6, family = plotFamily)
        plotRadial(relative, rownames(data), radial.lim = c(0, 0.3),
                   main = sprintf('Arginine codon usage across stages in %s\n', tissue))
        dev.off()
    }
}

plotRnaCodonDistribution <- function (codonUsageData, name, tissue) {
    on.exit(dev.off())
    pdf(sprintf('plots/usage/%s-codons-%s.pdf', name, tissue),
        width = 8, height = 3.5, family = plotFamily)

    for (aa in aminoAcids$Short) {
        codons <- rownames(subset(geneticCode, AA == aa))
        if (length(codons) == 1)
            next

        mrna <- codonUsageData[codons, grep(tissue, colnames(codonUsageData))]
        mrna <- relativeData(mrna[codons, ])

        long <- subset(aminoAcids, Short == aa)$Long

        colim <- sapply(stages, p(grep, colnames(mrna)))
        mrna <- `colnames<-`(mrna[, colim], stages)
        mrna <- mrna[, ncol(mrna) : 1]

        # We shuffle the order of the colours to avoid giving the misleading
        # impression of a gradient.
        grays <- gray.colors(length(colors) - 1)[c(3, 6, 2, 5, 1, 4, 7)]

        local({
            oldPar <- par(lty = 0, xpd = TRUE)
            on.exit(par(oldPar))
            barplot(mrna, horiz = TRUE, col = grays, axes = FALSE,
                    main = sprintf('Codon usage for %s', long),
                    xlab = 'Proportion of mRNA codon usage',
                    las = 1, names.arg = readable(colnames(mrna)))
            axis(1, 0 : 4 / 4, sprintf('%d%%', 0 : 4 * 25), cex.axis = 0.75)
            legendPos <- mrna[, ncol(mrna)]
            legendCol <- grays[legendPos != 0]
            legendPos <- cumsum(legendPos[legendPos != 0])
            legend <- names(legendPos)
            legendPos <- c(0, legendPos[-length(legendPos)])

            par(usr = c(0, 1, 0, 1))

            # Adjust label positions so they don’t overlap.
            legendWidths <- strwidth(legend)
            if (length(legend) > 1)
                for (i in 1 : (length(legend) - 1))
                    if (legendPos[i] + legendWidths[i] > legendPos[i + 1])
                        legendPos[i + 1] <- legendPos[i] + legendWidths[i]

            text(legendPos, 1, legend, pos = 4, col = legendCol,
                 cex = 0.75, offset = 0.1)
        })
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
    doPlot <- function (data, which) {
        map(.(data, name = {
            on.exit(dev.off())
            pdf(sprintf('plots/usage-sampling/codons-%s-%s.pdf', which, name))
            map(p(plotCodonUsageForAA, data, name), aminoAcids$Short)
        }), data, names(data))
    }

    map(doPlot, list(codonSampleMatrix, expressedCodonSampleMatrix),
        c('all', 'expressed')) %|% invisible
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
    map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/amino-acids-%s.pdf', name))
        plotAAUsage(data, readable(names(data)))
    }), list(codonSampleMatrix, expressedCodonSampleMatrix),
        c('all', 'expressed')) %|% invisible
}

if (! interactive()) {
    cat('# Generating mRNA codon usage data\n')
    mrnaLoadData()
    mrnaNormalizeData()
    generateCodonUsageData()
    loadAminoAcids()
    loadCodonMap()
    generateHighCodonUsageData()
    generateLowCodonUsageData()
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
    plotSpiderWeb('low-', lowCodonUsageData, lowAaUsageData)
    plotRnaCodonDistribution(stableCodonUsageData, 'stable', 'liver')
    plotRnaCodonDistribution(lowCodonUsageData, 'low', 'liver')

    cat('# Generate background mRNA codon usage plots\n')
    plotCodonBackground()

    cat('# Generate shuffled expression codon profiles\n')
    resampleCodonUsage()
    resampleExpressedCodonUsage()
    mkdir('plots/usage-sampling')
    plotCodonSampling()
    plotAminoAcidSampling()
}
