source('scripts/de.R')

countsForCondition <- function (conditions)
    trnaRawCounts[, colnames(trnaRawCounts) %in%
                    rownames(trnaMapping)[trnaMapping$Condition %in% conditions]]

runMeme <- function (foreground, background, outputPath) {
    memePath <- '/usr/local/Cellar/meme/4.9.0-p4/bin'
    memeBin <- file.path(memePath, 'meme')
    dustBin <- file.path(memePath, 'dust')
    filterFastaBin <- '../common/scripts/filter-fasta'
    markovModelBin <- '../common/scripts/markov-model-from-fasta'

    trnaUpstreamFastaFile <- '../common/data/trna-upstream-with-ids.fasta'

    mkdir(outputPath)

    system(sprintf('bash -c \'%s <(%s < %s %s) > %s\'',
                   dustBin,
                   filterFastaBin,
                   trnaUpstreamFastaFile,
                   paste(foreground, collapse = ' '),
                   file.path(outputPath, 'seq.fasta')))

    system(sprintf('bash -c \'%s <(%s <(%s < %s %s)) > %s\'',
                   markovModelBin,
                   dustBin,
                   filterFastaBin,
                   trnaUpstreamFastaFile,
                   paste(background, collapse = ' '),
                   file.path(outputPath, 'background.mm')))

    system(sprintf('%s %s -dna -oc %s -nostatus -maxsize 200000 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -revcomp -bfile %s',
                   memeBin,
                   file.path(outputPath, 'seq.fasta'),
                   file.path(outputPath, 'result'),
                   file.path(outputPath, 'background.mm')))
}

runMemeForStages <- function (tissue, a, b) {
    contrastStages <- paste(tissue, c(a, b), sep = '-')
    contrast <- paste(contrastStages, collapse = '/')
    deGenes <- rownames(trnaDeGenes[[contrast]])
    allGenes <- getExpressedtRNAs(countsForCondition(contrastStages))
    background <- setdiff(allGenes, deGenes)
    basePath <- file.path('results/meme', tissue, paste(a, b, sep = '-'))
    runMeme(deGenes, background, basePath)
}

runMemeOnAll <- function () {
    # Test expressed versus non expressed genes
    expressed <- rownames(trnaAnnotation)
    all <- rownames(trnaUnfilteredAnnotation)
    background <- setdiff(all, expressed)
    basePath <- 'results/meme/expressed'
    runMeme(expressed, background, basePath)
}

if (! interactive()) {
    trnaLoadData()
    trnaPairwiseDiffentialExpression()

    # For now, only look at e15.5/* contrasts

    for (tissue in tissues) {
        for (stage in stages[-1]) {
            runMemeForStages(tissue, 'e15.5', stage)
        }
    }

    runMemeOnAll()
}
