# Should be run as:
#
#   Rscript --vanilla common/scripts/generate-all.R
#   …

# Ensure we’re in the correct path …

argv <- commandArgs(trailingOnly = FALSE)
fopt <- '--file='
thisPath <- dirname(sub(fopt, '', argv[grep(fopt, argv)]))
setwd(file.path(thisPath, '..'))

runInLocal <- function (path, expr) {
    local({
        on.exit(setwd('../common'))
        setwd(file.path('..', path))
        evalq(expr)
    })
}

# RNA-seq data
runInLocal('rna', source('scripts/generate-all.R'))

# ChIP-seq data
runInLocal('chip', source('scripts/generate-all.R'))
