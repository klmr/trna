source('scripts/correlation.R')

plotSpiderWeb <- function () {
    plotRadial <- function (data, labels, main, ...) {
        layout(matrix(c(1, 2), 1, 2, byrow = TRUE), widths = c(1, 0.2))
        radial.plot(data, labels = labels, main = main,
                    line.col = c('#00000080', colors[indices(stages)]),
                    lwd = c(5, rep(2, ncol(data) - 1)), show.grid.labels = 3, ...)
        plot.new()
        legend('center', legend = readable(stages), fill = colors, border = NA,
               bty = 'n', xpd = NA)
    }

    # Enforce uniform order between tRNA and mRNA plots.
    isotypeData <- trnaByType[aminoAcids$Long, ]
    background <- table(trnaUnfilteredAnnotation$Type)
    background <- background[aminoAcids$Long]

    for (tissue in tissues) {
        data <- isotypeData[, grep(tissue, colnames(isotypeData))]
        relative <- relativeData(cbind(background, data))
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
    # Make row order consistent.
    rownames(data) <- revcomp(rownames(data))
    data <- data[codons, ]
    data[is.na(data)] <- 0

    main <- sprintf('Simulated %s isoacceptor abundance in %s',
                    long, readable(name))
    plotCodonBarplot(relativeData(data), main)
}

plotAcceptorSampling <- function () {
    doPlot <- function (data, which) {
        map(.(data, name = {
            on.exit(dev.off())
            pdf(sprintf('plots/usage-sampling/codons-%s-%s.pdf', which, name))
            map(p(plotAcceptorAbundanceForAA, data, name), aminoAcids$Short)
        }), data, names(data))
    }

    map(doPlot, list(acceptorSampleMatrix, expressedAcceptorSampleMatrix),
        c('all-genes', 'expressed-genes')) %|% invisible
}

plotIsotypeRadial <- function (data, col = colors[1 : ncol(data)], lwd = 2, main = as.character(substitute(data)))
    radial.plot(data, labels = rownames(data),
                line.col = col, lwd = lwd, show.grid.labels = 3,
                radial.lim = c(0, 0.1), main = main)

plotIsotypeUsage <- function (data, background, name) {
    n <- ncol(data[[1]])
    prepare <- cbind %|>% isotypeAbundance
    data <- do.call(cbind, map(prepare, data))
    data <- cbind(data, t(relativeData(matrix(background))))

    tcolors <- c(rep(transparent(colors[indices(stages)], 0.2), each = n),
                 colors['grey'])
    par(oma = rep(5, 4))
    plotIsotypeRadial(data, col = tcolors,
                      main = 'Isotype abundance with resampled expression')
}

plotIsotypeSampling <- function () {
    background <- table(trnaUnfilteredAnnotation$Type)
    background <- background[aminoAcids$Long]

    map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/amino-acids-%s.pdf', name))
        plotIsotypeUsage(data, background, readable(names(data)))
    }), list(acceptorSampleMatrix, expressedAcceptorSampleMatrix),
        c('all-genes', 'expressed-genes')) %|% invisible
}

plotRandomlyChosen <- function () {
    data <- pickRandomExpressions(acceptorSampleMatrix)
    plotIsotypeRadial(data, colors[indices(stages)], lwd = 2)
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

prepareSimulatedCodons <- function (codons, acceptors) {
    codons <- map(relativeData, codons)
    acceptors <- map(.(acc = {
        acc <- acceptorAbundance(acc)
        rownames(acc) <- revcomp(rownames(acc))

        # Ensure that non-used codons are present.
        codonsOnly <- setdiff(rownames(codons[[1]]), rownames(acc))
        codonNullRows <- do.call(rbind, map(.(. = rep(0, ncol(acc))), codonsOnly))
        acc <- rbind(acc, codonNullRows)
        # Make row order uniform with codons
        as.data.frame(acc[rownames(codons[[1]]), ])
    }), acceptors)

    list(codons = codons, acceptors = acceptors)
}

prepareSimulatedAas <- function (codons, acceptors) {
    prepareCodon <- function (data) {
        data <- groupby(data, geneticCode[rownames(data), 'AA'])
        # Enforce uniform oder between tRNA and mRNA plots.
        data <- data[aminoAcids$Short, ]
        rownames(data) <- aminoAcids[aminoAcids$Short == rownames(data), 'Long']
        relativeData(data)
    }

    prepareAcceptors <- function (data)
        isotypeAbundance(data)[aminoAcids$Long, ]

    list(codons = map(prepareCodon, codons),
         acceptors = map(prepareAcceptors, acceptors))
}

plotSimulatedDistribution <- function(type, codons, acceptors, datatype) {
    prepare <- list(`amino-acids` = prepareSimulatedAas,
                    codons = prepareSimulatedCodons)
    data <- prepare[[datatype]](codons, acceptors)

    # Calculate correlations for all tRNA simulations with the mean of the mRNA
    # codon simulations.

    codonMeans <- do.call(cbind, map(rowMeans, data$codons))
    conditions <- colnames(codonMeans)

    corr <- map(.(replicate = map(.(cond =
                cor(codonMeans[, cond], data$acceptors[[cond]][, replicate],
                    method = 'spearman')), conditions)),
                1 : ncol(data$acceptors[[1]]))

    pval <- map(.(replicate = map(.(cond =
                cor.test(codonMeans[, cond], data$acceptors[[cond]][, replicate],
                    method = 'spearman')$p.value), conditions)),
                1 : ncol(data$acceptors[[1]]))

    corr <- as.data.frame(do.call(rbind, map(unlist, corr)))
    pval <- as.data.frame(do.call(rbind, map(unlist, pval)))

    mkdir('results/usage-sampling')
    write.table(summary(corr), sep = '\t', quote = FALSE,
                sprintf('results/usage-sampling/%s-%s-correlations.tsv', datatype, type))

    write.table(summary(pval), sep = '\t', quote = FALSE,
                sprintf('results/usage-sampling/%s-%s-pvalues.tsv', datatype, type))

    # All boxes in one plot
    local({
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/%s-%s-correlations.pdf', datatype, type))
        boxplot(corr, pch = 16, las = 1,
                names = readable(colnames(corr)),
                main = sprintf('Correlation coefficients of simulated %s',
                               datatype))
    })

    # All data in one box
    local({
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/%s-%s-correlations-all-stages.pdf',
                    datatype, type))
        boxplot(unname(unlist(corr)), pch = 16, las = 1,
                main = sprintf('Correlation coefficients of simulated %s',
                               datatype))
    })

    # One plot per box, one box per stage
    map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage-sampling/%s-%s-correlation-%s.pdf', datatype, type, name))
        boxplot(data, pch = 16, las = 1,
                main = sprintf('Correlation coefficients of simulated %s in %s',
                               datatype, readable(name)))
    }), corr, colnames(corr)) %|% invisible
}

plotCodonAnticodonCorrelations <- function (trna, mrna, excludeZeros = FALSE,
                                            title, xlab, ylab) {
    maxima <- findPlotMaxima(trna, mrna)
    par(mfrow = c(4, 3))

    if (missing(title))
        title <- 'Codons in %s %s'

    if (missing(xlab))
        xlab <- 'Proportion of tRNA isoacceptors'

    if (missing(ylab))
        ylab <- 'Proportion of mRNA codon usage'

    map(.(tissue = {
        map(.(stage = {
            data <- columnsForCondition(trna, mrna, tissue, stage)
            colors <- if (excludeZeros)
                ifelse(data$trna == 0, last(colors), tissueColor[tissue]) else
                    tissueColor[tissue]
            plot(data$trna, data$mrna,
                 xlab = xlab, ylab = ylab,
                 main = sprintf(title, readable(stage), readable(tissue)),
                 xlim = c(0, maxima['trna']), ylim = c(0, maxima['mrna']),
                 col = colors, pch = 20, las = 1)
            if (excludeZeros)
                data <- data[data$trna != 0, ]
            model <- lm(mrna ~ trna, data)
            abline(model)
            par(usr = c(0, 1, 0, 1))
            rho <- cor(data$trna, data$mrna, method = 'spearman')
            r2 <- cor(data$trna, data$mrna, method = 'pearson')
            prho <- cor.test(data$trna, data$mrna, method = 'spearman')$p.value
            pr2 <- cor.test(data$trna, data$mrna, method = 'pearson')$p.value
            message(tissue, '-', stage, ': prho=', prho, ' pr2=', pr2)
            text(1, 0, bquote(atop(' ' ~ italic(p) == .(sprintf('%.2f', prho)) ~ (rho == .(sprintf('%.2f', rho))),
                                   italic(p) == .(sprintf('%.2f', pr2)) ~ (R^2 == .(sprintf('%.2f', r2))))),
                 adj = c(1.1, -0.1))

            rho
        }), stages) %|% unlist
    }), tissues)
}

plotCodonsByStage <- function (codonUsageData) {
    trnaCodons <- groupby(trnaNormDataCond, trnaAnnotation[rownames(trnaNormDataCond), 'Acceptor'])
    rownames(trnaCodons) <- revcomp(rownames(trnaCodons))
    onlym <- setdiff(rownames(codonUsageData), rownames(trnaCodons))
    trnaCodons[onlym, ] <- 0
    # Ensure same row order.
    trna <- relativeData(trnaCodons[rownames(codonUsageData), ])
    mrna <- relativeData(codonUsageData)
    plotCodonAnticodonCorrelations(trna, mrna, excludeZeros = TRUE)
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
        data <- map(.(tissue = {
            d <- data[, grep(tissue, colnames(data))]
            d <- d[, vapply(stages, p(grep, colnames(d)), numeric(1))]
            cor(d, method = 'spearman')
        }), tissues)
        data <- map(p(`diag<-`, NA), data)
        data <- map(p(`colnames<-`, stages), data)
        data <- map(p(`rownames<-`, stages), data)
        plotCountMatrix(data, title, '%0.2f', commonScale = globalScale)
    }

    #' @TODO Do not hard-code the value
    globalScale <- c(0.9, 1, 0.1)
    invisible(map(generateCorrelationPlot, names(methods), methods, titles))
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

        generateHighCodonUsageData()
        generateLowCodonUsageData()
        generateCodonBackgroundDist()
        generateCodonBackgroundUsage()
    })
    plotCodonsByType()

    corr <- map(.(data, name = {
        on.exit(dev.off())
        pdf(sprintf('plots/usage/%scodon-scatter.pdf', name),
            width = 7, height = 10, family = plotFamily)
        plotCodonsByStage(data)
    }), list(codonUsageData, stableCodonUsageData, lowCodonUsageData),
        c('', 'stable-', 'low-')) %|% p(unlist, recursive = FALSE)

    # Account for wobble positions in codon-anticodon pairing

    source('scripts/wobble-pairing-1.R')

    local({
        source('scripts/wobble-pairing-2.R')
        plotWobbleCorrelation()
    })

    # Count isoacceptor family sizes only

    local({
        source('scripts/family-sizes.R')
        plotIsoacceptorFamilySize()
        plotIsoacceptorFamilySizeNoWobbling()
    })

    local({
        on.exit(dev.off())
        pdf('plots/usage/codon-correlation-comparison.pdf')
        boxplot(corr, pch = 16, las = 1, border = tissueColor[names(corr)],
                names = rep(c('All', 'High', 'Low '), each = 2),
                main = 'Correlation coefficients between codon usage and isoacceptor abundance',
                ylim = c(0, 1))
    })

    local({
        on.exit(dev.off())
        pdf('plots/usage/amino-acid-scatter.pdf', width = 7, height = 10, family = plotFamily)
        plotAminAcidsByStage()
    })

    cat('# Generate shuffled expression isoacceptor profiles\n')
    resampleAcceptorAbundance()
    resampleExpressedAcceptorAbundance()
    mkdir('plots/usage-sampling')
    plotAcceptorSampling()
    plotIsotypeSampling()

    resampleCodonUsage()
    resampleExpressedCodonUsage()

    plotSimulatedDistribution('all-genes', codonSampleMatrix,
                              acceptorSampleMatrix, 'codons')
    plotSimulatedDistribution('expressed-genes', expressedCodonSampleMatrix,
                              expressedAcceptorSampleMatrix, 'codons')

    plotSimulatedDistribution('all-genes', codonSampleMatrix,
                              acceptorSampleMatrix, 'amino-acids')
    plotSimulatedDistribution('expressed-genes', expressedCodonSampleMatrix,
                              expressedAcceptorSampleMatrix, 'amino-acids')
}
