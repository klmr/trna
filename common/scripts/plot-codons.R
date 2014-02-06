plotCodonBarplot <- function (data, main) {
    # Copied from chip/scripts/correlation-plots.R:plotCodonsByType
    # Refactoring, yay.
    # We shuffle the order of the colours to avoid giving the misleading
    # impression of a gradient.
    grays <- gray.colors(length(colors) - 1)[c(3, 6, 2, 5, 1, 4, 7)]

    # Set xpd because the legend text at the top may overshoot plot area.
    oldPar <- par(xpd = TRUE)
    on.exit(par(oldPar))
    barplot(data, horiz = TRUE, col = grays, border = grays, axes = FALSE,
            xlab = 'Cumulative proportion of codons',
            ylab = 'Shuffled samples',
            space = 0, las = 1, names.arg = rep('', ncol(data)), main = main)
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
}
