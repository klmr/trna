radial.plot <- function (data, ...) {
    require(plotrix)
    plotrix::radial.plot(t(data), rp.type = '', ...)
    args <- list(t(data), rp.type = 'p', show.grid = FALSE,
                 show.radial.grid = FALSE, add = TRUE)
    args <- c(args, list(...))
    args$main <- NULL
    args$labels <- NULL
    do.call(plotrix::radial.plot, args)
}
