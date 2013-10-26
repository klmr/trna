source('scripts/correlation.R')

preparePdf <- function (name, ...)
    pdf(file.path('plots', 'correlation', paste(name, '.pdf', sep = '')), ...)

scatterPlot <- function (x, y, ...) {
    model <- lm(y ~ x)
    rho <- cor(x, y, method = 'spearman')
    points(x, y, ...)
    text(max(x), min(y),
         substitute(paste(rho, '=', r),
                    list(r = sprintf('%0.2f', rho))),
         adj = c(1, 0))
    #abline(0, 1, lty = 3)
    abline(model, lty = 2)
}

writeCodonUsageData <- function()
    write.table(codonUsageData, file = 'results/codon-usage.dat')

plotCorrelationMatrix <- function () {
    methods <- list(codon = codonUsageData, aa = aaUsageData)

    for (method in names(methods)) {
        data <- methods[[method]]
        preparePdf(sprintf('%s-usage', method), width = 8, height = 8)
        tissueScatterPlot <- function (t)
            function (x, y) {
                rows <- merged$Tissue == t
                scatterPlot(x[rows], y[rows],
                            cex = 0.8, pch = 16,
                            col = transparent(tissueColor[t]))
            }

        idx <- lapply(tissues, partial(grep, colnames(data)))
        only <- lapply(lapply(idx, function (i) data[, i]),
                       function(x) `colnames<-`(x, stages))
        # Create data set usable with `pairs` to plot different rows for top
        # and bottom half of matrix.
        merged <- do.call(rbind, only)
        merged$Tissue <- rep(tissues, each = nrow(data))

        pairs(merged[, -ncol(merged)],
              upper.panel = tissueScatterPlot('brain'),
              lower.panel = tissueScatterPlot('liver'))
        legend(0.38, 1.07 , legend = names(tissueColor), fill = tissueColor,
               border = NA, ncol = 2, xpd = NA, bty = 'n')
        dev.off()
    }
}

plotSpiderWeb <- function () {
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
        legend('center', legend = c(stages, 'back­\nground'), fill = tcolors, border = NA, bty = 'n', xpd = NA)
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
        pdf(sprintf('plots/usage/%s.pdf', tissue),
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
        pdf(sprintf('plots/usage/arginine-%s.pdf', tissue),
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

    cat('# Generating mRNA expression correlation plots\n')
    mkdir('plots/correlation')
    plotCorrelationMatrix()

    cat('# Generating mRNA codon usage plots\n')
    mkdir('plots/usage')
    plotSpiderWeb()
}
