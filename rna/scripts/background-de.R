source('../common/scripts/basic.R')
source('scripts/load-data.R')

# Calculate list of DE genes compared to background for each stage
# The background is formed by all other stages, respectively

mrnaDifferentFromBackground <- function (fdr = 0.05) {
    deGenesFromMean <- function (stage, tissue) {
        # Patch condition names in count data set
        condition <- mrnaMapping$Condition
        condB <- paste(tissue, stage, sep = '-')

        maskCondB <- boolmask(grep(condB, condition), length(condition))
        maskTissue <- boolmask(grep(tissue, condition), length(condition))
        condition[maskTissue & ! maskCondB] <- 'background'

        cds <- newCountDataSet(mrnaRawCounts, condition)
        cds <- estimateDispersions(estimateSizeFactors(cds))
        res <- nbinomTest(cds, 'background', condB)
        res[! is.na(res$padj) & res$padj <= fdr, ]
    }

    require(parallel)
    cores <- detectCores()
    sapply(tissues, lp(mclapply, stages, deGenesFromMean, mc.cores = cores),
                       simplify = FALSE)
}

reportDeGenesAcrossStages <- function () {
    de <- mrnaDifferentFromBackground()
    basepath <- 'results/de-background'

    outputDeSet <- function (stage, tissue) {
        filename = file.path(basepath, tissue, stage, ext = 'tsv')
        mkdir(file.path(basepath, tissue))
        write.table(de[[tissue]][[stage]], filename, quote = FALSE, sep = '\t',
                    col.names = NA)
    }

    invisible(sapply(tissues, lp(sapply, stages, outputDeSet)))
}

if (! interactive()) {
    mrnaLoadData()
    mrnaNormalizeData()

    reportDeGenesAcrossStages()
}
