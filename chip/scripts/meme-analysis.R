source('scripts/de.R')
source('scripts/expressed-per-stage.R')

memePath <- '/usr/local/Cellar/meme/4.9.0-p4/bin'
memeBin <- file.path(memePath, 'meme')
dustBin <- file.path(memePath, 'dust')
tomtomBin <- file.path(memePath, 'tomtom')
filterFastaBin <- '../common/scripts/filter-fasta'
markovModelBin <- '../common/scripts/markov-model-from-fasta'
trnaUpstreamFastaFile <- '../common/data/trna-upstream-with-ids.fasta'
memeDatabasePath <- '../common/data/meme/motif_databases'

meme <- function (foreground, background, outputPath) {
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
    extractHeader <- function (name, header)
        let(match = regexpr(sprintf('(?<=%s= )\\S+', name), header, perl = TRUE),
            as.numeric(regmatches(header, match)))

    file <- file.path(results, 'meme.txt')
    lines <- readLines(file)
    start <- grep('^\tMotif \\d+ position-specific probability matrix', lines) + 2
    length <- extractHeader('w', lines[start])
    nsites <- extractHeader('nsites', lines[start])
    evalue <- extractHeader('E', lines[start])
    matrix <- map(.(start, end = read.table(text = paste(lines[start : end],
                                                         collapse = '\n'),
                                            col.names = c('A', 'C', 'G', 'T'))),
                  start + 1, start + length)

    map(.(matrix, nsites, evalue =
          let(str = list(matrix = matrix, nsites = nsites, e = evalue),
              structure(str, class = 'pspm'))),
        matrix, nsites, evalue)
}

consensusSequence <- function (pspm)
    with(pspm, paste(colnames(matrix)[
                        vapply(as.data.frame(t(matrix)), which.max,
                               numeric(1))], collapse = ''))

pspmHeader <- function (pspm)
    paste(sprintf('MOTIF %s', consensusSequence(pspm)),
          sprintf('letter-probability matrix: alenght= %s w= %s nsites= %s E= %e',
                  ncol(pspm$matrix), nrow(pspm$matrix), pspm$nsites, pspm$e),
          sep = '\n\n')

print.pspm <- function (x, ...) {
    cat(pspmHeader(x), '\n')
    x
}

memeFileContent <- function (pspm, foreground) {
    # Format documented at
    # <http://meme.nbcr.net/meme/doc/meme-format.html#min_format>
    version <- 'MEME version 4.4\n'

    inputFasta <- local({
        on.exit(close(con))
        con <- pipe(sprintf('bash -c \'%s <(%s < %s %s)\'',
                            dustBin,
                            filterFastaBin,
                            trnaUpstreamFastaFile,
                            paste(foreground, collapse = ' ')))
        readLines(con)
    })

    # Remove Fasta headers and merge lines
    inputFasta <- paste(inputFasta[! grepl('^>', inputFasta)], collapse = '')
    # Normalise casing and remove Ns
    inputFasta <- gsub('N', '', toupper(inputFasta))
    freqs <- table(strsplit(inputFasta, ''))

    frequencies <- sprintf('Background letter frequencies:\n%s\n',
                           paste(mapply(paste, names(freqs), freqs), collapse = ' '))

    motif <- capture.output(write.table(pspm$matrix, col.names = FALSE,
                                        row.names = FALSE))

    c(version,
      # Alphabet omitted (can be inferred)
      # Strand omitted (can be inferred)
      frequencies,
      pspmHeader(pspm),
      motif)
}

tomtom <- function (pspm, foreground, outputPath, databases) {
    mkdir(outputPath)
    motifName <- consensusSequence(pspm)
    inputFile <- file.path(outputPath, motifName, ext = 'meme')
    writeLines(memeFileContent(pspm, foreground), inputFile)

    if (missing(databases))
        databases <- file.path(memeDatabasePath,
                               c('JASPAR_CORE_2009_vertebrates.meme',
                                 'uniprobe_mouse.meme'))

    system(sprintf('%s -no-ssc -oc %s -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10 %s %s',
                   tomtomBin,
                   file.path(outputPath, 'result', motifName),
                   inputFile,
                   databases))
}

getMemeSetupForStage <- function (tissue, a, b) {
    contrastStages <- paste(tissue, c(a, b), sep = '-')
    contrast <- paste(contrastStages, collapse = '/')
    deGenes <- rownames(trnaDeGenes[[contrast]])
    allGenes <- getExpressedtRNAs(countsForCondition(contrastStages))
    background <- setdiff(allGenes, deGenes)
    basePath <- file.path('results/meme', tissue, paste(a, b, sep = '-'))
    as.list(environment())
}

runMemeForStages <- function (tissue, a, b) {
    with(getMemeSetupForStage(tissue, a, b),
         meme(deGenes, background, basePath))
}

runTomtomForStage <- function (tissue, a, b) {
    with(getMemeSetupForStage(tissue, a, b), {
        pspm <- meme(deGenes, background, basePath)
        map(tomtom, pspm, list(deGenes), basePath) %|% invisible
    })
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
