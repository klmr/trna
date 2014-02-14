# Patch clipping region bug in plotrix::radial.plot

patchFunction  <- function () {
    funenv <- environment(plotrix::radial.plot)
    unlockBinding('radial.plot', funenv)

    unpatched <- body(plotrix::radial.plot)
    toPatch <- unpatched[[16]]
    toPatch[[2]] <- quote(NA)
    patched <- unpatched
    patched[[16]] <- toPatch
    f <- closure(formals(plotrix::radial.plot),
                 patched,
                 funenv)
    assign('radial.plot', f, envir = funenv)

    lockBinding('radial.plot', funenv)
}

patchFunction()

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
