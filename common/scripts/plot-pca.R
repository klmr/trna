# data = prcomp data
# dim1, dim2 = names of dimensions, e.g. 'Tissue', 'Stage'
# values1, values2 = values of dimensions
# mapping = data.frame with colnames at least c(dim1, dim2)
#
plotProgression <- function (data, dim1, dim2, values1, values2, mapping, main, cols, legend = TRUE) {
    symbols <- c(15, 16)
    # Distinct color for each stage
    makePalette <- function (cluster)
        colorRampPalette(c(lighten(cols[cluster], 0.7), cols[cluster]))(length(values2))
    palette <- list(makePalette(1), makePalette(2))
    names(palette) <- values1

    pchCluster <- function (ids)
        vapply(ids, function (id) symbols[grep(mapping[id, dim1], values1)], 0)

    # Note: `as.character` is required for factors, otherwise indices get confused.
    colStage <- function (ids)
        vapply(ids, function (id) palette[[as.character(mapping[id, dim1])]]
               [grep(mapping[id, dim2], values2)], '')

    variance <- summary(data)$importance['Proportion of Variance', ]
    axisLab <- '%s (%.0f%% variance explained)'

    lim <- sapply(1:2, function (p) c(min(data$rotation[, p]),
                                      max(data$rotation[, p])) * 1.5)

    par(bty = 'n', xpd = TRUE, las = 1)
    plot(data$rotation[, 1], data$rotation[, 2],
         xlab = sprintf(axisLab, 'PC1', variance[1] * 100),
         ylab = sprintf(axisLab, 'PC2', variance[2] * 100),
         col = colStage(rownames(data$rotation)),
         pch = pchCluster(rownames(data$rotation)),
         cex = 1.2, xlim = lim[, 1], ylim = lim[, 2],
         main = main)
    text(data$rotation[, 1], data$rotation[, 2],
         adj = c(-0.5, 0.5),
         labels = readable(as.character(mapping[rownames(data$rotation), dim2])))

    if (legend) {
        graphics::legend('topright', legend = readable(values1),
                         col = cols, inset = c(0, -max(data$rotation[, 2]) / 3),
                         pch = symbols, pt.lwd = 3, bty = 'n', cex = 1.2)
    }
}
