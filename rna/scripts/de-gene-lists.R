source('scripts/de.R')

mrnaWriteDeGenes <- function () {
    for (tissue in tissues) {
        base <- file.path('results/de', tissue)
        mkdir(base)

        contrasts <- names(mrnaDeGenes)
        contrasts <- contrasts[grepl(tissue, contrasts)]

        map(.(contrast =
              let(filename = file.path(base, gsub(paste0(tissue, '-'), '',
                        sub('/', '-', contrast)), ext = 'tsv'),
                  write.table(mrnaDeGenes[[contrast]], filename,
                              col.names = NA, sep = '\t', quote = FALSE))),
            contrasts)
    }
}

if (! interactive()) {
    cat('# Generating protein-coding DE gene lists\n')
    mrnaLoadData()
    mrnaPairwiseDifferentialExpression()
    mrnaWriteDeGenes()
}
