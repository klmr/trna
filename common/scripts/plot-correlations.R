#' Plot pairwise correlations for a matrix of data points
plotCorrelationsFor <- function (data, main = as.character(substitute(data))) {
    plotter <- function (tissue)
        function (i, j, ...) {
            cols <- sapply(c(i, j), p(grep, colnames(data)) %.%
                           lp(paste, tissue, sep = '-'))
            plotSampleIdx <- if (useSample)
                sample.int(nrow(data), 2000) else 1 : nrow(data)
            plotSample <- data[plotSampleIdx, cols]

            plot(plotSample, col = col[tissue], pch = 16, cex = 0.5,
                 log = if (useSample) 'xy' else '', axes = FALSE,
                 bty = 'n', xlab = '', ylab = '', xlim = lablim, ylim = lablim)

            if (j == stages[1]) axis(3)
            if (j == last(stages)) axis(1)
            if (i == stages[1]) axis(2, las = 1)
            if (i == last(stages)) axis(4, las = 1)

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
    col <- (if (nrow(data) > 200) {
            opacity <- max(200 / nrow(data), 0.05)
            p(transparent, opacity)
        } else id)(tissueColor)
    names(col) <- names(tissueColor)

    # Use uniform axis limits.
    lablim <- c(0, max(unlist(data)))

    axislabel <- sprintf(if (useSample) 'log(%s)' else '%s', 'Gene expression')
    plotPairwise(stages, plotter('liver'), plotter('brain'), diagonal,
                 main = main, xlab = axislabel, ylab = axislabel)
}
