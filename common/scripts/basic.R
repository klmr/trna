# Silence package loading – ensure it’s only run once:

if (! exists('oldLibrary')) {
    oldLibrary <- library
    library <- function (package)
        suppressPackageStartupMessages(
            oldLibrary(as.character(substitute(package)),
                       verbose = FALSE, character.only = TRUE))

    oldRequire <- require
    require <- function (package)
        suppressPackageStartupMessages(
            oldRequire(as.character(substitute(package)), character.only = TRUE))
}

# Qualified path to make this usable everywhere.
source('../common/scripts/definitions.R', chdir = TRUE)
