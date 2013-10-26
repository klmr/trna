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
                          ...) {
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
            upper(axes[i], axes[j], ...)
        else if (i > j)
            lower(axes[i], axes[j], ...)
        else
            do.call(diagonal, c(list(axes[i]), .diag))
    }

    .par$mar <- .par$mar %else% rep(1, 4)
    .par$mfrow <- rep(length(axes), 2)
    oldPar <- do.call(par, .par)
    on.exit(par(oldPar))

    dots <- list(...)
    for (name in names(dots))
        .diag[[name]] <- .diag[[name]] %else% dots[[name]]

    idx <- indices(axes)
    apply(expand.grid(idx, idx), ROWS, callPlot)
}

preparePdf <- function (name, ...)
    pdf(file.path('plots', 'correlation', name, ext = '.pdf'), ...)

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
