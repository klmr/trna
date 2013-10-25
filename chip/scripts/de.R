source('../common/scripts/basic.R')
source('scripts/load-data.R')

trnaPairwiseDiffentialExpression <- function () {
    if (exists('trnaDeGenes'))
        return()

    threshold <- 0.05

    cond <- function (x) sprintf('%s-%s', tissue, stages[x])
    contrast <- function (i, j) sprintf('%s/%s', cond(i), cond(j))

    # Declares stageList[stage][otherStage] = NULL for all stages
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
                result <- nbinomTest(trnaCds, cond(i), cond(j))
                significant <- subset(result, ! is.na(padj) & padj < threshold)
                deResults[[tissue]][[stages[i]]][[stages[j]]] <- result
                deGenes[[contrast(i, j)]] <- significant
                # The matrix is symmetric
                deCounts[[tissue]][i, j] <- nrow(significant)
                deCounts[[tissue]][j, i] <- nrow(significant)
            }
        }
    }

    trnaDeResults <<- deResults
    trnaDeGenes <<- deGenes
    trnaDeCounts <<- deCounts
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
