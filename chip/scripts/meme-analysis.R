source('scripts/de.R')
source('scripts/expressed-per-stage.R')

meme <- function (foreground, background, outputPath) {
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

    system(sprintf('%s %s -dna -oc %s -maxsize 200000 -mod zoops -nmotifs 3 -minw 6 -maxw 50 -revcomp -bfile %s',
                   memeBin,
                   file.path(outputPath, 'seq.fasta'),
                   file.path(outputPath, 'result'),
                   file.path(outputPath, 'background.mm')))

    parsePspm(file.path(outputPath, 'result'))
}

parsePspm <- function (results) {
    # Parse output raw text file and retrieve PSSMs for all motifs.
    # Note: we could also parse the XML file using XPath but the XML output file
    # is quite frankly not very nice. Parsing the raw text is easier. Fail.
    file <- file.path(results, 'meme.txt')
    lines <- readLines(file)
    start <- grep('^\tMotif \\d+ position-specific probability matrix', lines) + 3
    header <- start - 1
    lengthMatch <- regexpr('(?<=w= )\\d+', lines[header], perl = TRUE)
    length <- as.numeric(regmatches(lines[header], lengthMatch))
    map(.(start, end = read.table(text = paste(lines[start : end], collapse = '\n'),
                                  col.names = c('A', 'C', 'G', 'T'))),
        start, start + length - 1)
}

runMemeForStages <- function (tissue, a, b) {
    contrastStages <- paste(tissue, c(a, b), sep = '-')
    contrast <- paste(contrastStages, collapse = '/')
    deGenes <- rownames(trnaDeGenes[[contrast]])
    allGenes <- getExpressedtRNAs(countsForCondition(contrastStages))
    background <- setdiff(allGenes, deGenes)
    basePath <- file.path('results/meme', tissue, paste(a, b, sep = '-'))
    meme(deGenes, background, basePath)
}

runMemeOnAll <- function () {
    # Test expressed versus non expressed genes
    expressed <- rownames(trnaAnnotation)
    all <- rownames(trnaUnfilteredAnnotation)
    background <- setdiff(all, expressed)
    basePath <- 'results/meme/expressed'
    meme(expressed, background, basePath)
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
