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

generateStableCodonUsageData <- function (low, high) {
    getCodonUsageAndReportProcess <- function (c, method) {
        range <- quantile(mrnaNormDataCond[, c], c(low, high))
        data <- mrnaNormDataCond[, c, drop = FALSE]
        data <- data[data >= range[1] & data <= range[2], , drop = FALSE]
        on.exit(cat('.'))
        codonUsage(method, data)
    }

    conditions <- unique(mrnaMapping$Condition)

    getCodonDataFrame <- function (method)
        as.data.frame(lapply(conditions, getCodonUsageAndReportProcess,
                             method = method))

    usageData <- map(getCodonDataFrame, c('codon', 'aa'))
    map(p(`colnames<-`, conditions), usageData)
}

generateHighCodonUsageData <- function () {
    # “Stable” codon usage reflects the codon usage not of all genes, but
    # instead only considers stably expressed genes in a given tissue, i.e.
    # those expressed within a certain range.
    # This is an attempt to smooth out high variance predominantly observed in
    # lowly expressed genes.
    # Conversely, we also exclude the very highly expressed genes in order to
    # de-emphasise extreme outliers in the data. To achieve this, we use the
    # data between the 90th and 95th quantile.

    cat('Generating stable codon usage data')
    usageData <- generateStableCodonUsageData(0.9, 0.95)
    cat('done\n')

    stableCodonUsageData <<- usageData[[1]]
    stableAaUsageData <<- usageData[[2]]
}

generateLowCodonUsageData <- function () {
    # Same as high stable usage, but for lowly expressed genes (but excluding
    # the lowest, highly variable ones) as a control.

    cat('Generating stable codon usage data')
    usageData <- generateStableCodonUsageData(0.25, 0.5)
    cat('done\n')

    lowCodonUsageData <<- usageData[[1]]
    lowAaUsageData <<- usageData[[2]]
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

#'@TODO Use `mrnaDoResample` to remove code redundancy
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

resampleExpressedCodonUsage <- function () {
    if (exists('expressedCodonSampleMatrix'))
        return()

    sampledCodonsFile <- '../common/cache/sampled-expressed-codons.RData'
    mkdir('../common/cache')
    suppressWarnings(loaded <- tryCatch(load(sampledCodonsFile), error = .(e = '')))

    if (! isTRUE(all.equal(loaded, 'expressedCodonSampleMatrix'))) {
        cat('Generating codon usage for permuted expressions')

        expressed <- getExpressedmRNAs(mrnaRawCounts)
        # Unlike for tRNA, we do not re-scale the library size now.
        # The reason is that we simply want a filtered set of expressed mRNA
        # genes as comparison to the total set, having the same expression values.
        mrnaNDC <- mrnaNormDataCond[expressed, ]
        expressedCodonSampleMatrix <- mrnaDoResample(mrnaNDC, 'liver')
        save(expressedCodonSampleMatrix, file = sampledCodonsFile)
    }

    expressedCodonSampleMatrix <<- expressedCodonSampleMatrix
}

getExpressedmRNAs <- function (counts) {
    # Filter out mRNAs which are unexpressed across all conditions.
    # We call "unexpressed" any mRNA whose expression value is below
    # a set threshold theta for at least one replicate.
    # Return a vector of the names of expressed mRNAs

    dos <- function (cond)
        rownames(subset(mrnaMapping, Condition == cond))
    mrnas <- function (cond)
        apply(meetsThreshold[, dos(cond)], ROWS, all)

    # Value determined by examining the data. Yes, it also works for mRNA.
    threshold <- 10
    conditions <- unique(mrnaMapping[colnames(counts), 'Condition'])

    meetsThreshold <- counts > threshold
    meetsThresholdPerCond <- sapply(conditions, mrnas)
    expressed <- apply(meetsThresholdPerCond, ROWS, any)
    names(expressed[expressed])
}

mrnaDoResample <- function (data, tissue, samples = 100) {
    data <- data[, grep(tissue, colnames(data))]

    require(parallel)

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
    codonSampleMatrix
}
