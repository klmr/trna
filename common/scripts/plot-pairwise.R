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
plotPairwise <- function (axes, upper, lower, diagonal,
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

    title(main, sub, xlab, ylab)
}
