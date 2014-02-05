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
    read.table(text = resultData, row.names = 1, header = FALSE)
}

loadCodonMap <- function ()
    geneticCode <<- read.table('../common/data/genetic_code.tsv', header = FALSE,
                               row.names = 1, col.names = c('', 'AA'),
                               stringsAsFactors = FALSE)

generateCodonUsageData <- function () {
    if (exists('codonUsageData'))
        return()

    codonUsageFile <- '../common/cache/codons-usage.RData'
    mkdir('../common/cache')

    suppressWarnings(loaded <- tryCatch(load(codonUsageFile), error = .(e = '')))

    if (! isTRUE(all.equal(loaded, 'usageData'))) {
        cat('Cache file not found -- re-generating codon usage data.')

        getCodonUsageAndReportProcess <- function (c, method) {
            on.exit(cat('.'))
            codonUsage(method, mrnaNormDataCond[, c, drop = FALSE])
        }

        getCodonDataFrame <- function (method)
            as.data.frame(lapply(conditions, getCodonUsageAndReportProcess,
                                 method = method))

        conditions <- unique(mrnaMapping$Condition)
        usageData <- lapply(c('codon', 'aa'), getCodonDataFrame)
        for (i in indices(usageData))
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
    # The background distribution assumes uniform expression strength for all
    # genes, so we just set all expression values to 1.
    if (exists('codonBackgroundDist'))
        return()

    allGenes <- data.frame(rep(1, nrow(mrnaNormDataCond)),
                           row.names = rownames(mrnaNormDataCond))
    codonBackgroundDist <<- codonUsage('codon', allGenes)
    aaBackgroundDist <<- codonUsage('aa', allGenes)
}

backgroundCodonUsage <- function (strand, frame) {
    script <- 'scripts/codon_background_usage.py'
    transcripts <- '../common/data/Mus_musculus.NCBIM37.67.transcripts.fa'

    resultData <- paste(readcmd('python', script, transcripts, strand, frame),
                        collapse = '\n')
    read.table(text = resultData, row.names = 1, col.names = c('', 'Count'))
}

generateCodonBackgroundUsage <- function () {
    if (exists('overallCodonBackground'))
        return()
    strand <- c(1, -1)
    frameshift <- 0 : 2
    overallCodonBackground <-
        apply(expand.grid(strand, frameshift), ROWS,
              .(row = backgroundCodonUsage(row[1], row[2])))
    overallAaBackground <-
        map(.(x = groupby(x, geneticCode[rownames(x), 1])),
            overallCodonBackground)
    overallCodonBackground <<- do.call(cbind, overallCodonBackground)
    overallAaBackground <<- map(.(x = x[rownames(x) != 'Stop', , drop = FALSE]),
                                overallAaBackground) %|% lp(do.call, cbind)
}

shuffleRows <- function (df)
    `rownames<-`(df[sample.int(nrow(df)), ], rownames(df))

resampleCodonUsage <- function () {
    if (exists('codonSampleMatrix'))
        return()

    sampledCodonsFile <- '../common/cache/sampled-codons.RData'
    mkdir('../common/cache')
    suppressWarnings(loaded <- tryCatch(load(sampledCodonsFile), error = .(e = '')))

    if (! isTRUE(all.equal(loaded, 'codonSampleMatrix'))) {
        samples <- 100

        data <- mrnaNormDataCond[, grep('liver', colnames(mrnaNormDataCond))]

        require(parallel)

        cat('Generating codon usage for permuted expressions')
        codonSamples <- mclapply(1 : samples, .(i = {
            d <- shuffleRows(data)
            on.exit(cat('.'))
            do.call(cbind, sapply(colnames(d),
                                  lp(`[`, d) %|>% lp(codonUsage, 'multi_codon')))
        }), mc.cores = detectCores())
        cat('\n')

        # Extract all same conditions
        codonSampleMatrix <- map(.(col = map(.(m = m[, col]), codonSamples) %|%
                                   lp(do.call, cbind)), 1 : ncol(data))

        codonSampleMatrix <- map(p(`rownames<-`, rownames(codonUsageData)),
                                 codonSampleMatrix)
        names(codonSampleMatrix) <- colnames(data)
        save(codonSampleMatrix, file = sampledCodonsFile)
    }

    codonSampleMatrix <<- codonSampleMatrix
}
