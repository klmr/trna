# Silence package loading

library <- function (package)
    suppressPackageStartupMessages(
        base::library(as.character(substitute(package)),
                      verbose = FALSE, character.only = TRUE))

require <- function (package)
    suppressPackageStartupMessages(
       base::require(as.character(substitute(package)), character.only = TRUE))

# Qualified path to make this usable everywhere.
source('../common/scripts/definitions.R', chdir = TRUE)
