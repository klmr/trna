plotCountMatrix <- function (counts, main) {
    maxval <- lapply(counts, partial(max, na.rm = TRUE))
    hsvliver <- rgb2hsv(col2rgb(gray.colors(maxval$liver + 1, start = 0.95, end = 0.3)))
    hsvbrain <- rgb2hsv(col2rgb(gray.colors(maxval$brain + 1, start = 0.95, end = 0.3)))

    # We want white--red so we ramp saturation with the inverse of value and leave hue (0 = red).
    hsvliver['s', ] <- 1 - hsvliver['v', ]
    hsvliver['v', ] <- 1
    hsvbrain['s', ] <- 1 - hsvbrain['v', ]
    hsvbrain['h', ] <- 0.13 # Orange-brown-yellowish
    hsvbrain['v', ] <- 0.95
    colors <- list(liver = hsv2col(hsvliver), brain = hsv2col(hsvbrain))

    w <- nrow(counts[[1]])
    par(mfrow = c(w, w), mar = rep(0, 4), oma = c(1, 1, 4, 10))
    for (i in seq(w)) {
        for (j in seq(w)) {
            tissue <- if (i < j) 'liver' else 'brain'
            value <- counts[[tissue]][i, j]
            plot.new()
            rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4],
                 col = colors[[tissue]][value + 1], border = NA)
            if (i == j)
                text(0.5, 0.5, labels = colnames(counts[[tissue]])[i], cex = 2)
            else
                text(0.5, 0.5, labels = value, cex = 1.5)
        }
    }

    title(main, cex.main = 2, outer = TRUE)

    # Legend
    op <- par(usr = c(-1.1, -0.1, -4.9, -3.9), xpd = NA)
    xpos = list(liver = 0.25, brain = 0.75)
    for (tissue in names(xpos)) {
        x <- xpos[[tissue]]
        points(rep(x, maxval[[tissue]]), (1 : maxval[[tissue]]) / maxval[[tissue]],
               col = colors[[tissue]], pch = 15, cex = 3)
        text(c(x, x), c(0, 1), labels = c(0, maxval[[tissue]]), pos = 4, offset = 1)
        text(x, 0, labels = sprintf('  %s', capitalize(tissue)), pos = 1, offset = 1)
    }
    par(op)
}
