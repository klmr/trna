# Dataset specific settings hard-coded here {{{
trnaDataFile <- '../common/data/trna-merged.tsv'
trnaMappingFile <- '../common/data/chip-sample-map.tsv'
trnaAnnotationFile <- '../common/data/tRNA.tsv'
# }}}

# Data setup {{{
trnaLoadData <- function () {
    if (exists('trnaRawCounts'))
        return()

    trnaRawCounts <- read.table(trnaDataFile, header = TRUE, row.names = 1)

    trnaMapping <- read.table(trnaMappingFile, header = FALSE, row.names = 1)
    colnames(trnaMapping) <- c('AB', 'Tissue', 'Stage')
    rownames(trnaMapping) <- sprintf('do%s', rownames(trnaMapping))
    trnaMapping <- trnaMapping[trnaMapping$AB == 'PolIII', ]
    trnaMapping$AB <- NULL
    trnaMapping$Condition <- paste(trnaMapping$Tissue, trnaMapping$Stage, sep = '-')
    # Ensure corresponding order.
    trnaMapping <- trnaMapping[colnames(trnaRawCounts), ]

    trnaAnnotation <- read.table(trnaAnnotationFile, header = FALSE)
    rownames(trnaAnnotation) <- sprintf('%s%s', trnaAnnotation[[1]], strstrip(trnaAnnotation[[2]]))
    trnaAnnotation <- trnaAnnotation[, c(3, 6, 7, 4, 5)]
    colnames(trnaAnnotation) <- c('Iso', 'Type', 'Acceptor', 'Start', 'End')
    trnaAnnotation$Strand <- ifelse(trnaAnnotation$Start < trnaAnnotation$End, '+', '-')
    trnaAnnotation$Iso <- strstrip(trnaAnnotation$Iso)
    trnaAnnotation$Type <- strstrip(trnaAnnotation$Type)
    trnaAnnotation$Chr <- sub('chr', '',
                              regmatches(rownames(trnaAnnotation),
                                         regexpr('chr[^.]+', rownames(trnaAnnotation))))

    # Swap inverted start and end …

    for (i in 1 : nrow(trnaAnnotation)) {
        if (trnaAnnotation[i, 'Strand'] == '-') {
            tmp <- trnaAnnotation[i, 'Start']
            trnaAnnotation[i, 'Start'] <- trnaAnnotation[i, 'End']
            trnaAnnotation[i, 'End'] <- tmp
        }
    }

    # Ensure corresponding order.
    trnaRawCounts <- trnaRawCounts[rownames(trnaAnnotation), ]

    trnaMapping <<- trnaMapping

    # Filter out unexpressed tRNAs
    require(DESeq) # Necessary for `counts` even though it’s loaded anyway.
    normCounts <- counts(trnaGetCountDataSet(trnaRawCounts), normalized = TRUE)
    expressed <- getExpressedtRNAs(normCounts)

    trnaUnfilteredRawCounts <<- trnaRawCounts
    trnaUnfilteredAnnotation <<- trnaAnnotation
    trnaRawCounts <<- trnaRawCounts[expressed, ]
    trnaAnnotation <<- trnaAnnotation[expressed, ]
}

getExpressedtRNAs <- function (counts) {
    # Filter out tRNAs which are unexpressed across all conditions.
    # We call "unexpressed" any tRNA whose expression value is below
    # a set threshold theta for at least one replicate.
    # Return a vector of the names of expressed tRNAs

    dosForCondition <- function (cond)
        rownames(subset(trnaMapping, Condition == cond))
    tRNAsInCondition <- function (cond)
        apply(meetsThreshold[, dosForCondition(cond)], ROWS, all)

    # Value determined by examining the data.
    conditions <- unique(trnaMapping[colnames(counts), 'Condition'])

    meetsThreshold <- counts > trnaExpressedThreshold
    meetsThresholdPerCond <- sapply(conditions, tRNAsInCondition)
    trnasExpressed <- apply(meetsThresholdPerCond, ROWS, any)
    names(trnasExpressed[trnasExpressed])
}
# }}}

# Normalized count data {{{
trnaGetCountDataSet <- function (data) {
    require(DESeq)
    cds <- newCountDataSet(data, trnaMapping[colnames(data), 'Condition'])
    estimateSizeFactors(cds)
}

trnaSetupCountDataSet <- function () {
    if (exists('trnaCds'))
        return()

    trnaCds <<- estimateDispersions(trnaGetCountDataSet(trnaRawCounts))
}

trnaMergeReplicates <- function (normData) {
    conditions <- unique(trnaMapping$Condition)

    merged <- matrix(nrow = nrow(normData), ncol = length(conditions))
    colnames(merged) <- conditions
    rownames(merged) <- rownames(normData)

    for (cond in conditions) {
        replicates = normData[, rownames(trnaMapping[trnaMapping$Condition == cond, ])]
        if (! is.data.frame(replicates))
            merged[, cond] <- replicates
        else
            merged[, cond] <- rowMeans(replicates)
    }

    as.data.frame(merged)
}

trnaNormalizeData <- function () {
    if (exists('trnaNormData'))
        return()

    trnaNormData <<- as.data.frame(counts(trnaCds, normalized = TRUE))
    trnaNormDataCond <<- trnaMergeReplicates(trnaNormData)
}
# }}}

mkdir('results')
mkdir('plots')
