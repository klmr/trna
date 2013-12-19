source('../common/scripts/basic.R')
source('scripts/load-data.R')

pairwiseDifferentialExpression <- function (data, threshold, dispEst = estimateDispersions) {
    cond <- function (x) sprintf('%s-%s', tissue, stages[x])
    contrast <- function (i, j) sprintf('%s/%s', cond(i), cond(j))

    stageList <- map(.(map(.(NULL), stages)), stages)
    deResults <- list(liver = stageList, brain = stageList)
    deGenes <- list()

    # Declares a pairwise matrix of counts for all stages
    sigC <- structure(rep(NA, length(stages) ^ 2),
                      .Dim = c(length(stages), length(stages)),
                      .Dimnames = list(stages, stages))
    deCounts <- list(liver = sigC, brain = sigC)

    for (tissue in tissues) {
        for (i in 1 : (length(stages) - 1)) {
            for (j in (i + 1) : length(stages)) {
                contrastData <- data[, rownames(subset(trnaMapping,
                                                       Condition %in% c(cond(i), cond(j))))]
                cds <- dispEst(estimateSizeFactors(trnaGetCountDataSet(contrastData)))
                result <- nbinomTest(cds, cond(i), cond(j))
                significant <- subset(result, ! is.na(padj) & padj < threshold)
                deResults[[tissue]][[stages[i]]][[stages[j]]] <- result
                deGenes[[contrast(i, j)]] <- significant
                # The matrix is symmetric
                deCounts[[tissue]][i, j] <- nrow(significant)
                deCounts[[tissue]][j, i] <- nrow(significant)
            }
        }
    }

    list(results = deResults, de = deGenes, counts = deCounts)
}

trnaPairwiseDiffentialExpression <- function () {
    if (exists('trnaDeGenes'))
        return()

    results <- pairwiseDifferentialExpression(trnaRawCounts, 0.05)

    trnaDeResults <<- results$results
    trnaDeGenes <<- results$de
    trnaDeCounts <<- results$counts
}

trnaDetailedDe <- function () {
    if (exists('trnaAccDe'))
        return()

    dispEst <- c(Acceptor = estimateDispersions,
                 Type = fun(cds = {
                            on.exit(assign('lp', lpartial, globalenv()))
                            rm(lp, envir = environment(lp))
                            estimateDispersions(cds, fitType = 'local')
                        }))

    getData <- function (x) groupby(trnaRawCounts, trnaAnnotation[[x]])

    data <- sapply(c('Acceptor', 'Type'), getData)
    trnaAccDe <<- pairwiseDifferentialExpression(data$Acceptor, 0.05)
    trnaTypeDe <<- pairwiseDifferentialExpression(data$Type, 0.05)
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
