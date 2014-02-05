source('scripts/correlation.R')

plotSpiderWeb <- function () {
    plotRadial <- function (data, labels, main, ...) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 0.2))
        radial.plot(data, labels = labels, main = main, line.col = colors,
                    lwd = 2, show.grid.labels = 3, ...)
        plot.new()
        legend('center', legend = readable(stages), fill = colors, border = NA,
               bty = 'n', xpd = NA)
    }

    # Enforce uniform order between tRNA and mRNA plots.
    isotypeData <- trnaByType[aminoAcids$Long, ]

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

plotAcceptorAbundanceForAA <- function (aa, data, name) {
    codons <- rownames(subset(geneticCode, AA == aa))
    if (length(codons) == 1)
        return()

    long <- subset(aminoAcids, Short == aa)$Long
    isotypes <- rownames(subset(trnaAnnotation, Type == long))
    data <- groupby(data[isotypes, ], trnaAnnotation[isotypes, 'Acceptor'])
    rownames(data) <- revcomp(rownames(data))

    main <- sprintf('Simulated %s isoacceptor abundance in %s',
                    long, readable(name))
    plotCodonBarplot(relativeData(data), main)
}

plotAcceptorSampling <- function () {
    map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/codons-%s.pdf', name))
        map(p(plotAcceptorAbundanceForAA, data, name), aminoAcids$Short)
    }), acceptorSampleMatrix, names(acceptorSampleMatrix)) %|% invisible
}

plotIsotypeUsage <- function (data, name) {
    plotRadial <- function (data, labels, main, ...) {
        tcolors <- c(rep('#40404050', ncol(data) - 1), colors[1])
        lwd <- c(rep(2, ncol(data) - 1), 5)
        radial.plot(data, labels = labels, main = main,
                    line.col = tcolors, lwd = lwd, show.grid.labels = 3, ...)
    }

    data <- groupby(data, trnaAnnotation$Type)
    # Enforce uniform oder between tRNA and mRNA plots.
    data <- data[aminoAcids$Long, ]
    plotRadial(relativeData(data), rownames(data),
               radial.lim = c(0, 0.1),
               main = sprintf('Isotype abundance with resampled expression for %s\n', name))
}

plotIsotypeSampling <- function () {
    allBackground <- trnaNormDataCond[, grep('liver', colnames(trnaNormDataCond))]
    allBackground <- allBackground[, vapply(stages, p(grep, colnames(allBackground)), numeric(1))]
    map(.(data, background, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/amino-acids-%s.pdf', name))
        plotIsotypeUsage(cbind(data, background), readable(name))
    }), acceptorSampleMatrix, allBackground, names(acceptorSampleMatrix)) %|% invisible
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

            # Reverse the order of the stages.
            mrna <- mrna[, ncol(mrna) : 1]
            trna <- trna[, ncol(trna) : 1]

            # We shuffle the order of the colours to avoid giving the misleading
            # impression of a gradient.
            grays <- gray.colors(length(colors) - 1)[c(3, 6, 2, 5, 1, 4, 7)]

            local({
                oldPar <- par(mfrow = c(1, 2), lty = 0, xpd = TRUE)
                on.exit(par(oldPar))
                mapply(function (data, title) {
                        barplot(data, horiz = TRUE, col = grays, axes = FALSE,
                                xlab = sprintf('Proportion of %s', title),
                                las = 1, names.arg = readable(colnames(data)))
                        axis(1, 0 : 4 / 4, sprintf('%d%%', 0 : 4 * 25), cex.axis = 0.75)
                        legendPos <- data[, ncol(data)]
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
            prho <- cor.test(data$trna, data$mrna, method = 'spearman')$p.value
            pr2 <- cor.test(data$trna, data$mrna, method = 'pearson')$p.value
            message(tissue, '-', stage, ': prho=', prho, ' pr2=', pr2)
            text(1, 0, bquote(atop(' ' ~ italic(p) == .(sprintf('%.2f', prho)) ~ (rho == .(sprintf('%.2f', rho))),
                                   italic(p) == .(sprintf('%.2f', pr2)) ~ (R^2 == .(sprintf('%.2f', r2))))),
                 adj = c(1.1, -0.1))
        }
    }
}

plotChipCorrelation <- function () {
    methods <- list(gene = trnaNormDataCond,
                    acceptor = trnaByAcceptor,
                    type = trnaByType)
    titles <- list(gene = 'tRNA gene expression correlation',
                   acceptor = 'tRNA isoacceptor family correlation',
                   type = 'tRNA isotype correlation')
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

    mrnaLoadData()
    mrnaNormalizeData()
    local({
        oldwd <- getwd(); on.exit(setwd(oldwd))
        setwd('../rna')
        generateCodonUsageData()
    })
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

    cat('# Generate shuffled expression isoacceptor profiles\n')
    resampleAcceptorAbundance()
    mkdir('plots/usage-sampling')
    plotAcceptorSampling()
    plotIsotypeSampling()
}
