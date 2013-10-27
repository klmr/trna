source('scripts/correlation.R')

#' Prepares a pairs plot and invokes a callback for each actual plot
#'
#' This function differs from \link{\code{pairs}} in that it only sets up the
#' plot window, it doesn’t itself manage the data. This makes it more flexible.
#' @param axes vector of labels denoting combinations to plot
#' @param upper function that is called for plots above the diagonal
#' @param lower function that is called for plots below the diagonal
#'      (defaults to \code{upper})
#' @param diagonal function that is called for plots on the diagonal
#'      (defaults to a function which just prints the labels)
#' @param .par list of additional \link{\code{par}}ameters to set
#' @param .diag additional parameters to pass to \code{diagonal}
#' @param main see \link{\code{title}}
#' @param sub see \link{\code{title}}
#' @param xlab see \link{\code{title}}
#' @param ylab see \link{\code{title}}
#' @param ... remaining arguments are passed to each plot function
#' @note The \code{upper} and \code{lower} callback function receive the
#'  (\code{i}, \code{j}) entries of the \code{axes} corresponding to its plot.
#'  The \code{diagonal} function receives its corresponding \code{axes} label.
#' @examples
#'  plotPairs(1:5, function (i, j) plot(x, y[i, j]))
plotPairwise <- function (axes,
                          upper,
                          lower,
                          diagonal,
                          .par = list(),
                          .diag = list(cex = 3),
                          main = as.character(substitute(axes)),
                          sub = '', xlab = 'x', ylab = 'y', ...) {
    if (missing(lower))
        lower <- upper
    if (missing(diagonal))
        diagonal <- function (x, ...) {
            image(matrix(1), col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n')
            text(0, 0, x, ...)
        }


    callPlot <- function (ij) {
        i <- ij[1]
        j <- ij[2]
        if (i < j)
            lower(axes[i], axes[j], ...)
        else if (i > j)
            upper(axes[i], axes[j], ...)
        else
            do.call(diagonal, c(list(axes[i]), .diag))
    }

    dots <- list(...)
    for (name in names(dots))
        .diag[[name]] <- .diag[[name]] %else% dots[[name]]

    .par$mar <- .par$mar %else% rep(1, 4)
    .par$oma <- .par$oma %else% c(5, 5, 5, 0)
    .par$mfrow <- rep(length(axes), 2)

    local({
        oldPar <- do.call(par, .par)
        on.exit(par(oldPar))
        idx <- indices(axes)
        apply(expand.grid(idx, idx), ROWS, callPlot)
    })

    title(main, '', xlab, ylab)
}

plotCorrelationsFor <- function (data, main = as.character(substitute(data))) {
    plotter <- function (tissue)
        function (i, j, ...) {
            cols <- sapply(c(i, j), p(grep, colnames(data)) %.%
                           lp(paste, tissue, sep = '-'))
            plotSampleIdx <- if (useSample)
                sample.int(nrow(data), 2000) else 1 : nrow(data)
            plotSample <- data[plotSampleIdx, cols]

            col <- (if (useSample) p(transparent, 0.05) else id)(tissueColor[tissue])
            plot(plotSample, col = col, pch = 20, cex = 0.5,
                 log = if (useSample) 'xy' else '',
                 bty = 'n', xlab = '', ylab = '')

            par(xlog = FALSE, ylog = FALSE, usr = c(0, 1, 0, 1))
            rho <- cor(data[, cols[1]], data[, cols[2]], method = 'spearman')
            text(0.5, 0.5, sprintf('%0.2f', rho), cex = 1.5)
        }

    diagonal <- function (stage, ...) {
        image(matrix(1), col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n')
        text(0, 0, readable(stage), cex = 2, font = 2)
    }

    #' Sample if > 2000 data points.
    useSample <- nrow(data) > 2000

    axislabel <- sprintf(if (useSample) 'log(%s)' else '%s', 'Gene expression')
    plotPairwise(stages, plotter('liver'), plotter('brain'), diagonal,
                 main = main, xlab = axislabel, ylab = axislabel)
}

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

    cat('# Generating mRNA expression & usage correlation plots\n')
    mkdir('plots/correlation')
    plotRnaCorrelation()

    cat('# Generating mRNA codon usage plots\n')
    mkdir('plots/usage')
    plotSpiderWeb()
}
