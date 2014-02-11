source('../common/scripts/basic.R')
source('scripts/load-data.R')
local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/usage.R')
})

trnaGroupFamilyAndType <- function () {
    trnaByAcceptor <<- groupby(trnaNormDataCond, trnaAnnotation$Acceptor, sum)
    trnaByType <<- groupby(trnaNormDataCond, trnaAnnotation$Type, sum)
}

resampleAcceptorAbundance <- function () {
    samples <- 100
    cds <- trnaGetCountDataSet(trnaUnfilteredRawCounts)
    data <- trnaMergeReplicates(as.data.frame(counts(cds, normalized = TRUE)))
    data <- data[, grep('liver', colnames(data))]
    # Reorder columns by stages since they are unordered
    data <- data[, vapply(stages, grep, numeric(1), colnames(data))]

    sampleColumn <- function (col)
        map(.(. = shuffleRows(data[, col, drop = FALSE])), 1 : samples) %|%
        lp(do.call, cbind)
    acceptorSampleMatrix <- setNames(map(sampleColumn, colnames(data)),
                                     colnames(data))

    acceptorSampleMatrix <<- acceptorSampleMatrix
}

acceptorAbundance <- function (trnaAbundance) {
    annotation <-
        if(nrow(trnaAbundance) == nrow(trnaAnnotation)) trnaAnnotation else
            trnaUnfilteredAnnotation

    relativeData(groupby(trnaAbundance, annotation$Acceptor))
}

isotypeAbundance <- function (trnaAbundance) {
    annotation <-
        if(nrow(trnaAbundance) == nrow(trnaAnnotation)) trnaAnnotation else
            trnaUnfilteredAnnotation

    data <- groupby(trnaAbundance, annotation$Type)
    # Enforce uniform order between tRNA and mRNA plots.
    data <- data[aminoAcids$Long, ]
    relativeData(data)
}

pickRandomExpressions <- function (acceptorSampleMatrix) {
    # From the sampled tRNA gene expression values, pick random values for each
    # amino acid to simulate an expression profile. Then check whether it
    # corresponds closely to the actually observed amino acid usage.

    pickRandom <- function (data) {
        aas <- sample(1 : ncol(data), nrow(data))
        idx <- cbind(indices(aas), aas)
        distribution <- apply(idx, ROWS, .(i = data[i[1], i[2]]))
        setNames(distribution / sum(distribution), rownames(data))
    }

    # Do it for each condition separately.
    do.call(cbind, map(isotypeAbundance %|>% pickRandom, acceptorSampleMatrix))
}
