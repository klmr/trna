plotCountMatrix <- function (counts, main, format = '%s', useSameScale = FALSE) {
    cols <- 100
    gradient <- function (tissue)
        colorRampPalette(c('white', tissueColor[tissue]),
                         interpolate = 'spline')(cols + 1)
    colors <- map(gradient, tissues)

    range <- if (useSameScale) {
        r <- range(counts, na.rm = TRUE)
        setNames(rep(list(c(r, r[2] - r[1])), 2), names(counts))
    } else {
        r <- map(p(range, na.rm = TRUE), counts)
        map(.(x = c(x, x[2] - x[1])), r)
    }

    w <- nrow(counts[[1]])
    par(mfrow = c(w, w), mar = rep(0, 4), oma = c(1, 1, 4, 10))
    for (i in seq(w)) {
        for (j in seq(w)) {
            tissue <- if (i < j) 'liver' else 'brain'
            value <- counts[[tissue]][i, j]
            color <- (value - range[[tissue]][1]) / range[[tissue]][3] * cols
            plot.new()
            rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4],
                 col = colors[[tissue]][color], border = NA)
            if (i == j)
                text(0.5, 0.5,
                     labels = readable(colnames(counts[[tissue]])[i]),
                     cex = 2, font = 2)
            else
                text(0.5, 0.5, labels = sprintf(format, value), cex = 2)
        }
    }

    title(main, cex.main = 2, outer = TRUE)

    # Legend
    op <- par(usr = c(-1.1, -0.1, -4.9, -3.9), xpd = NA)
    xpos = list(liver = 0.25, brain = 0.75)
    for (tissue in names(xpos)) {
        x <- xpos[[tissue]]
        points(rep(x, cols), (1 : cols) / cols,
               col = colors[[tissue]], pch = 15, cex = 3)
        text(c(x, x), c(0, 1), labels = sprintf(format, range[[tissue]][-3]), pos = 4, offset = 1)
        text(x, 0, labels = sprintf('  %s', capitalize(tissue)), pos = 1, offset = 1)
    }
    par(op)
}
