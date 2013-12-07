source('../common/scripts/basic.R')
source('scripts/load-data.R')
source('scripts/download-transcripts.R')

codonUsage <- function (method, countData) {
    script <- sprintf('scripts/%s_frequency.py', method)
    # If the following file doesn't exist, run
    #
    #   Rscript scripts/download-transcripts.R
    #
    transcripts <- '../common/data/Mus_musculus.NCBIM37.67.transcripts.fa'
    if (! file.exists(transcripts)) {
        warning('Transcript file not found. Downloading gene transcripts. This takes a VERY long time.',
                call. = FALSE, immediate. = TRUE)
        downloadTranscripts(transcripts)
        cat('done\n')
    }

    tmp <- tempfile()
    on.exit(unlink(tmp))
    write.table(countData, file = tmp, sep = '\t', col.names = FALSE, quote = FALSE)
    resultData <- paste(readcmd('python', script, transcripts, '<', tmp),
                        collapse = '\n')
    read.table(text = resultData, row.names = 1, col.names = c('', 'Count'))
}

loadCodonMap <- function ()
    geneticCode <<- read.table('../common/data/genetic_code.tsv', header = FALSE,
                               row.names = 1, col.names = c('', 'AA'))

generateCodonUsageData <- function () {
    if (exists('codonUsageData'))
        return()

    codonUsageFile <- '../common/cache/codons-usage.RData'
    mkdir('../common/cache')

    suppressWarnings(loaded <- tryCatch(load(codonUsageFile), error = .('')))

    if (! isTRUE(all.equal(loaded, 'usageData'))) {
        cat('Cache file not found -- re-generating codon usage data.')

        getCodonUsageAndReportProcess <- function (c, method) {
            result <- codonUsage(method, mrnaNormDataCond[, c, drop = FALSE])
            cat('.')
            result
        }

        getCodonDataFrame <- function (method)
            as.data.frame(lapply(conditions, getCodonUsageAndReportProcess,
                                 method = method))

        conditions <- unique(mrnaMapping$Condition)
        usageData <- lapply(c('codon', 'aa'), getCodonDataFrame)
        for (i in 1 : length(usageData))
            colnames(usageData[[i]]) <- conditions

        cat('done\n')
        save(usageData, file = codonUsageFile)
    }

    codonUsageData <<- usageData[[1]]
    aaUsageData <<- usageData[[2]]
}

generateStableCodonUsageData <- function () {
    # “Stable” codon usage reflects the codon usage not of all genes, but
    # instead only considers highly expressed genes in a given tissue.
    # This is an attempt to smooth out high variance predominantly observed in
    # lowly expressed genes.
    # Conversely, we also exclude the very highly expressed genes in order to
    # de-emphasise extreme outliers in the data. To achieve this, we use the
    # data between the 90th and 95th quantile.

    conditions <- unique(mrnaMapping$Condition)

    cat('Generating stable codon usage data for.')

    getCodonUsageAndReportProcess <- function (c, method) {
        range <- quantile(mrnaNormDataCond[, c], c(0.9, 0.95))
        data <- mrnaNormDataCond[, c, drop = FALSE]
        data <- data[data >= range[1] & data <= range[2], , drop = FALSE]
        on.exit(cat('.'))
        codonUsage(method, data)
    }

    getCodonDataFrame <- function (method)
        as.data.frame(lapply(conditions, getCodonUsageAndReportProcess,
                             method = method))

    usageData <- lapply(c('codon', 'aa'), getCodonDataFrame)
    for (i in indices(usageData))
        colnames(usageData[[i]]) <- conditions

    cat('done\n')

    stableCodonUsageData <<- usageData[[1]]
    stableAaUsageData <<- usageData[[2]]
}

generateCodonBackgroundDist <- function () {
    if (exists('codonBackgroundDist'))
        return()

    allGenes <- data.frame(rep(1, nrow(mrnaNormDataCond)),
                           row.names = rownames(mrnaNormDataCond))
    codonBackgroundDist <<- codonUsage('codon', allGenes)
    aaBackgroundDist <<- codonUsage('aa', allGenes)
}
