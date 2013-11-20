source('../common/scripts/basic.R')

plotColorScheme <- function () {
    on.exit(dev.off())
    pdf('../common/color-scheme.pdf', width = 8, height = 4)
    par(mfrow = c(2, 1), mar = c(3, 1, 3, 1))
    image(as.matrix(indices(colors)), col = colors, xaxt = 'n', yaxt = 'n',
          bty = 'n', main = 'tRNA paper colour scheme')
    par(usr = c(0, 1, 0, 1))
    pos <- indices(colors) / length(colors)
    text(pos - pos[1] / 2, 0.75, colors, adj = c(0.5, 0.5))
    text(pos[3 : 4] - pos[1] / 2, 0.25, readable(names(tissueColor)),
         adj = c(0.5, 0.5))

    image(as.matrix(indices(contrastColors)), col = contrastColors,
          xaxt = 'n', yaxt = 'n', bty = 'n', main = 'Contrast colours')
} 

if (! interactive()) {
    plotColorScheme()
}
