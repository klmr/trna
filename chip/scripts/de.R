source('../common/scripts/basic.R')
source('scripts/load-data.R')

pairwiseDifferentialExpression <- function (data, mapping, threshold) {
    require(DESeq2)

    cond <- function (x) sprintf('%s-%s', tissue, stages[x])
    contrast <- function (i, j) sprintf('%s/%s', cond(i), cond(j))

    stageList <- map(.(map(.(NULL), stages)), stages)
    results <- list(liver = stageList, brain = stageList)
    de <- list()

    # Declares a pairwise matrix of counts for all stages
    sigC <- structure(rep(NA, length(stages) ^ 2),
                      .Dim = c(length(stages), length(stages)),
                      .Dimnames = list(stages, stages))
    counts <- list(liver = sigC, brain = sigC)

    maxProgress <- length(tissues) * (length(stages) ^ 2 - length(stages)) / 2
    currentProgress <- 0
    for (tissue in tissues) {
        for (i in 1 : (length(stages) - 1)) {
            for (j in (i + 1) : length(stages)) {
                progress(currentProgress, maxProgress)
                currentProgress <- currentProgress + 1
                testContrast <- c(cond(i), cond(j))
                testLibraries <- rownames(subset(mapping,
                                                 Condition %in% testContrast))
                contrastData <- data[, testLibraries]
                condition <- factor(mapping[testLibraries, 'Condition'],
                                    testContrast)
                cds <- DESeqDataSetFromMatrix(contrastData,
                                              data.frame(condition = condition),
                                              ~ condition)
                cds <- suppressMessages(DESeq(cds))
                result <- results(cds)
                significant <- subset(result, ! is.na(padj) & padj < threshold)
                results[[tissue]][[stages[i]]][[stages[j]]] <- result
                de[[contrast(i, j)]] <- significant
                # The matrix is symmetric
                counts[[tissue]][i, j] <- nrow(significant)
                counts[[tissue]][j, i] <- nrow(significant)
            }
        }
    }
    progress(currentProgress, maxProgress)

    list(results = results, de = de, counts = counts)
}

trnaPairwiseDiffentialExpression <- function () {
    if (exists('trnaDeGenes'))
        return()

    results <- pairwiseDifferentialExpression(trnaRawCounts, trnaMapping, 0.05)

    trnaDeResults <<- results$results
    trnaDeGenes <<- results$de
    trnaDeCounts <<- results$counts
}

trnaDetailedDe <- function () {
    if (exists('trnaAccDe'))
        return()

    getData <- function (x) groupby(trnaRawCounts, trnaAnnotation[[x]])

    data <- map(getData, c('Acceptor', 'Type'))
    trnaAccDe <<- pairwiseDifferentialExpression(data$Acceptor, trnaMapping, 0.05)
    trnaTypeDe <<- pairwiseDifferentialExpression(data$Type, trnaMapping, 0.05)
}

trnaTissueDifferentialExpression <- function () {
    if (exists('trnaTissueDeGenes'))
        return()

    doDe <- function (stage) {
        cond <- p(paste, stage, sep = '-')
        nbinomTest(trnaCds, cond('liver'), cond('brain'))
    }
    trnaTissueDeResults <<- lapply(stages, p(sorted, 'padj') %.% doDe)
}

deGenesForContrast <- function (contrast, threshold = 0.05)
    contrast[contrast$padj < threshold & ! is.na(contrast$padj), ]

nonDeGenesForContrast <- function (contrast, threshold = 0.05)
    contrast[setdiff(rownames(contrast),
                     rownames(deGenesForContrast(contrast, threshold))), ]
