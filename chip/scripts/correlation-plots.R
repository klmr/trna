source('scripts/correlation.R')

# TODO Redundant definition, consolidate in common helper file.
preparePdf <- function (name, ...)
    pdf(file.path('plots/correlation', name, ext = '.pdf'), family = plotFamily, ...)

pointSpread <- function (x, y, ...) {
    points(x, y, ...)
    linm <- lm(y ~ x)
    rho <- cor(x, y, method = 'spearman')
    text(max(x), min(y), bquote(rho == .(sprintf('%.2f', rho))), adj = c(1, 0))
    #abline(0, 1, lty = 3)
    abline(linm, lty = 2)
}

ps <- partial(pointSpread, cex = 0.8, pch = 16, col = transparent(colors[1]))

# Override radial.plot options

radial.plot <- function (data, ...) {
    require(plotrix)
    plotrix::radial.plot(t(data), rp.type = '', ...)
    args <- list(t(data), rp.type = 'p', show.grid = FALSE,
                 show.radial.grid = FALSE, add = TRUE)
    args <- c(args, list(...))
    args$main <- NULL
    args$labels <- NULL
    do.call(plotrix::radial.plot, args)
}

plotSpiderWeb <- function () {
    plotRadial <- function (data, labels, main, ...) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 0.2))
        radial.plot(data, labels = labels, main = main, line.col = colors,
                    lwd = 2, show.grid.labels = 3, ...)
        plot.new()
        legend('center', legend = readable(stages), fill = colors, border = NA,
               bty = 'n', xpd = NA)
    }

    isotypeData <- groupby(trnaNormDataCond, trnaAnnotation$Type)
    # Enforce uniform order between tRNA and mRNA plots.
    isotypeData <- isotypeData[aminoAcids$Long, ]

    for (tissue in tissues) {
        data <- isotypeData[, grep(tissue, colnames(isotypeData))]
        relative <- relativeData(data)
        pdf(sprintf('plots/usage/trna-%s.pdf', tissue),
            width = 6, height = 6, family = plotFamily)
        plotRadial(relative, rownames(data),
                   radial.lim = c(0, 0.1),
                   main = sprintf('Isotype occupancy across stages in %s\n', tissue))
        dev.off()
    }
}

plotCodonsByType <- function () {
    for (tissue in tissues) {
        pdf(sprintf('plots/usage/codons-%s.pdf', tissue),
            width = 8, height = 3.5, family = plotFamily)

        for (aa in aminoAcids$Short) {
            codons <- rownames(subset(geneticCode, AA == aa))
            if (length(codons) == 1)
                next

            mrna <- codonUsageData[codons, grep(tissue, colnames(codonUsageData))]

            long <- subset(aminoAcids, Short == aa)$Long
            isotypes <- rownames(subset(trnaAnnotation, Type == long))
            trna <- groupby(trnaNormDataCond[isotypes, grep(tissue, colnames(trnaNormDataCond))],
                            trnaAnnotation[isotypes, 'Acceptor'])
            rownames(trna) <- revcomp(rownames(trna))

            cat(sprintf('%s (%s)\n', long, aa))
            cat(codons)
            cat('\n\n')

            # Make them have common row names.

            onlym <- setdiff(rownames(mrna), rownames(trna))
            onlyt <- setdiff(rownames(trna), rownames(mrna))

            if (length(onlyt) > 0)
                mrna[onlyt, ] <- 0
            if (length(onlym) > 0)
                trna[onlym, ] <- 0

            # Make sure their rows are in the same order, and make relative

            mrna <- relativeData(mrna[codons, ])
            trna <- relativeData(trna[codons, ])

            # Finally, make sure the columns are in the same order.
            colim <- sapply(stages, p(grep, colnames(mrna)))
            colit <- sapply(stages, p(grep, colnames(trna)))

            mrna <- `colnames<-`(mrna[, colim], stages)
            trna <- `colnames<-`(trna[, colit], stages)

            local({
                oldPar <- par(mfrow = c(1,2), lty = 0, xpd = TRUE)
                on.exit(par(oldPar))
                mapply(function (data, title) {
                        barplot(data, horiz = TRUE, col = colors, axes = FALSE,
                                xlab = sprintf('Proportion of %s', title),
                                las = 1, names.arg = readable(stages))
                        axis(1, 0 : 4 / 4, sprintf('%d%%', 0 : 4 * 25), cex.axis = 0.75)
                        legendPos <- data[, ncol(data)]
                        legendCol <- colors[legendPos != 0]
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
                    },
                    list(trna, mrna),
                    c('tRNA isoacceptor occupancy', 'mRNA codon usage'))
            })
            title(main = sprintf('Codon usage for %s', long), xpd = NA)
        }

        dev.off()
    }
}

# Scatter plot of amino acids by stages

columnsForCondition <- function (trna, mrna, tissue, stage) {
    cols <- lapply(list(trna, mrna), colnames)
    lc <- length(cols[[1]])
    thistissue <- lapply(cols, p(boolmask, lc) %.% lp(grep, tissue))
    thisstage <- lapply(cols, p(boolmask, lc) %.% lp(grep, stage))
    cols <- mapply(`&`, thisstage, thistissue)
    data.frame(trna = trna[, cols[, 1]],
               mrna = mrna[, cols[, 2]])
}

findPlotMaxima <- function (trna, mrna) {
    maxima <- list()

    # Determine maxima to unify all axis limits. Messy
    for (tissue in tissues)
        for (stage in stages)
            maxima[[length(maxima) + 1]] <-
                apply(columnsForCondition(trna, mrna, tissue, stage), COLS, max)

    apply(do.call(rbind, maxima), COLS, max) * 1.05
}

plotCodonsByStage <- function () {
    trnaCodons <- groupby(trnaNormDataCond, trnaAnnotation[rownames(trnaNormDataCond), 'Acceptor'])
    rownames(trnaCodons) <- revcomp(rownames(trnaCodons))
    onlym <- setdiff(rownames(codonUsageData), rownames(trnaCodons))
    trnaCodons[onlym, ] <- 0
    # Ensure same row order.
    trna <- relativeData(trnaCodons[rownames(codonUsageData), ])
    mrna <- relativeData(codonUsageData)
    maxima <- findPlotMaxima(trna, mrna)

    par(mfrow = c(4, 3))

    for (tissue in tissues) {
        for (stage in stages) {
            data <- columnsForCondition(trna, mrna, tissue, stage)
            plot(data$trna, data$mrna,
                 xlab = 'Proportion of tRNA isoacceptors',
                 ylab = 'Proportion of mRNA codon usage',
                 main = sprintf('Codons in %s %s', readable(stage), readable(tissue)),
                 xlim = c(0, maxima['trna']), ylim = c(0, maxima['mrna']),
                 col = ifelse(data$trna == 0, last(colors), tissueColor[tissue]),
                 pch = 20, las = 1)
            cd <- data[data$trna != 0, ]
            model <- lm(mrna ~ trna, cd)
            abline(model)
            par(usr = c(0, 1, 0, 1))
            rho <- cor(cd$trna, cd$mrna, method = 'spearman')
            r2 <- cor(cd$trna, cd$mrna, method = 'pearson')
            text(1, 0, bquote(atop(' ' ~ rho == .(sprintf('%.2f', rho)),
                                      R^2 == .(sprintf('%.2f', r2)))),
                 adj = c(1.1, -0.1))
        }
    }
}

plotAminAcidsByStage <- function () {
    trnaTypeUsage <- groupby(trnaNormDataCond, trnaAnnotation[rownames(trnaNormDataCond), 'Type'])
    mrnaTypeUsage <- groupby(codonUsageData, geneticCode[rownames(codonUsageData), 'AA'])
    mrnaTypeUsage <- mrnaTypeUsage[! rownames(mrnaTypeUsage) == 'Stop', ]
    rownames(mrnaTypeUsage) <- sapply(rownames(mrnaTypeUsage),
                                      function (aa) subset(aminoAcids, Short == aa)$Long)

    # Ensure same row order.
    trna <- relativeData(trnaTypeUsage[aminoAcids$Long, ])
    mrna <- relativeData(mrnaTypeUsage[aminoAcids$Long, ])
    maxima <- findPlotMaxima(trna, mrna)

    par(mfrow = c(4, 3))

    for (tissue in tissues) {
        for (stage in stages) {
            data <- columnsForCondition(trna, mrna, tissue, stage)
            plot(data$trna, data$mrna,
                 xlab = 'Proportion of tRNA isotypes',
                 ylab = 'Proportion of mRNA amino acid usage',
                 main = sprintf('Amino acids in %s %s', readable(stage), readable(tissue)),
                 xlim = c(0, maxima['trna']), ylim = c(0, maxima['mrna']),
                 pch = 20, las = 1, col = tissueColor[tissue])
            model <- lm(mrna ~ trna, data)
            ci <- predict(model, interval = 'confidence', level = 0.99)
            ciHigh <- sort(ci[, 'upr'])
            ciLow <- sort(ci[, 'lwr'])
            ordx <- sort(data$trna)
            xx <- c(ordx, rev(ordx), ordx[1])
            yy <- c(ciLow, rev(ciHigh), ciLow[1])
            abline(model)
            outliers <- (data$mrna < ci[, 'lwr']) | (data$mrna > ci[, 'upr'])
            text(data$trna, data$mrna, labels = ifelse(outliers, rownames(data), ''), pos = 4)
            polygon(xx, yy, col = '#00000040', border = NA)
            par(usr = c(0, 1, 0, 1))
            rho <- cor(data$trna, data$mrna, method = 'spearman')
            r2 <- cor(data$trna, data$mrna, method = 'pearson')
            text(1, 0, bquote(atop(' ' ~ rho == .(sprintf('%.2f', rho)),
                                   R^2 == .(sprintf('%.2f', r2)))),
                 adj = c(1.1, -0.1))
        }
    }
}

plotChipCorrelation <- function () {
    methods <- list(gene = trnaNormDataCond,
                    acceptor = trnaByAcceptor)
    titles <- list(gene = 'tRNA gene expression correlation',
                   acceptor = 'tRNA isoacceptor family correlation')
    generateCorrelationPlot <- function (name, data, title) {
        on.exit(dev.off())
        pdf(file.path('plots', 'correlation',
                      sprintf('%s-usage', name), ext = 'pdf'))
        plotCorrelationsFor(data, title)
    }

    invisible(mapply(generateCorrelationPlot, names(methods), methods, titles))
}

if (! interactive()) {
    cat('# Generating correlation plots\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()
    trnaGroupFamilyAndType()

    mkdir('plots/correlation')
    plotChipCorrelation()

    loadAminoAcids()
    loadGeneticCode()
    mkdir('plots/usage')
    plotSpiderWeb()

    generateCodonUsageData()
    plotCodonsByType()
    local({
        on.exit(dev.off())
        pdf('plots/usage/codon-scatter.pdf', width = 7, height = 10, family = plotFamily)
        plotCodonsByStage()
    })
    local({
        on.exit(dev.off())
        pdf('plots/usage/amino-acid-scatter.pdf', width = 7, height = 10, family = plotFamily)
        plotAminAcidsByStage()
    })
}
