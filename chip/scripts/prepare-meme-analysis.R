source('scripts/de.R')

countsForCondition <- function (conditions)
    trnaRawCounts[, colnames(trnaRawCounts) %in%
                    rownames(trnaMapping)[trnaMapping$Condition %in% conditions]]

if (! interactive()) {
    trnaLoadData()
    trnaPairwiseDiffentialExpression()

    filterFasta <- '../common/scripts/filter-fasta'
    markovModel <- '../common/scripts/markov-model-from-fasta'

    trnaUpstreamFastaFile <- '../common/data/tRNA-upstream-with-ids.fasta'

    # For now, only look at e15.5/* contrasts

    for (tissue in tissues) {
        for (stage in stages[-1]) {
            contrastStages <- paste('liver', c('e15.5', stage), sep = '-')
            contrast <- paste(contrastStages, collapse = '/')
            deGenes <- rownames(trnaDeGenes[[contrast]])
            allGenes <- getExpressedtRNAs(countsForCondition(contrastStages))
            background <- setdiff(allGenes, deGenes)

            basePath <- file.path('results/meme', tissue, stage)
            mkdir(basePath)

            system(sprintf('%s < %s %s > %s',
                           filterFasta,
                           trnaUpstreamFastaFile,
                           paste(deGenes, collapse = ' '),
                           file.path(basePath, 'seq.fasta')))

            system(sprintf('bash -c \'%s <(%s < %s %s) > %s\'',
                           markovModel,
                           filterFasta,
                           trnaUpstreamFastaFile,
                           paste(background, collapse = ' '),
                           file.path(basePath, 'background.mm')))
        }
    }
}
