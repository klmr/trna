source('../common/scripts/basic.R')
source('scripts/load-data.R')

mrnaPairwiseDifferentialExpression <- function () {
    if (exists('mrnaDeGenes'))
        return()

    deGenesFile <- '../common/cache/de-genes.RData'
    mkdir('../common/cache')

    suppressWarnings(loaded <- tryCatch(load(deGenesFile), error = .('')))

    if (! isTRUE(all.equal(loaded, c('mrnaDeGenes', 'mrnaDeResults')))) {
        cat('Cache file not found -- re-generating DE gene lists.')

        contrastFor <- function (t, s) sprintf('%s-%s', t, s)
        threshold <- 0.05
        # Declares stageList[stage][otherStage] = NULL for all stages
        stageList <- map(.(map(.(NULL), stages)), stages)
        mrnaDeGenes <- list(liver = stageList, brain = stageList)
        mrnaDeResults <- mrnaDeGenes

        for (tissue in tissues) {
            for (i in 1 : (length(stages) - 1)) {
                for (j in (i + 1) : length(stages)) {
                    a <- contrastFor(tissue, stages[i])
                    b <- contrastFor(tissue, stages[j])
                    result <- nbinomTest(mrnaCds, a, b)
                    de <- subset(result, ! is.na(padj) & padj < threshold)
                    mrnaDeResults[[tissue]][[stages[i]]][[stages[j]]] <- result
                    mrnaDeGenes[[tissue]][[stages[i]]][[stages[j]]] <- de
                    cat('.')
                }
            }
        }

        cat('done\n')
        save(mrnaDeGenes, mrnaDeResults, file = deGenesFile)
    }

    sigC <- matrix(nrow = length(stages), ncol = length(stages))
    mrnaDeCount <- list(liver = sigC, brain = sigC)

    for (tissue in tissues) {
        for (i in 1 : (length(stages) - 1)) {
            for (j in (i + 1) : length(stages)) {
                de <- mrnaDeGenes[[tissue]][[stages[i]]][[stages[j]]]
                mrnaDeCount[[tissue]][i, j] <- nrow(de)
                mrnaDeCount[[tissue]][j, i] <- nrow(de)
            }
        }
    }

    setNames <- function (cnt) `colnames<-`(`rownames<-`(cnt, stages), stages)
    mrnaDeCount <- lapply(mrnaDeCount, setNames)

    mrnaDeResults <<- mrnaDeResults
    mrnaDeGenes <<- mrnaDeGenes
    mrnaDeCount <<- mrnaDeCount
}
