source('../common/scripts/basic.R')
source('scripts/load-data.R')
source('scripts/de.R')

local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/de.R')
})

require(plyr) # for adply

inNeighborhood <- function (gene, enhancers, windowSize)
    count(enhancers$Chr == gene$Chr &
          enhancers$Pos >= gene$Start - windowSize &
          enhancers$Pos <= gene$End + windowSize)

countEnhancersPerGene <- function (genes, enhancers, windowSize)
    adply(genes, ROWS, inNeighborhood, enhancers, windowSize)$V1

deGenes <- function (allGenes, threshold)
    subset(allGenes, ! is.na(padj) & padj <= threshold)

genePositions <- function (genes, annotation)
    annotation[rownames(genes), c('Chr', 'Start', 'End')]

activePositions <- function (condition) {
    activeGeneNames <- getExpressedtRNAs(countsForCondition(condition))
    activeGenes <- trnaUnfilteredRawCounts[activeGeneNames, ]
    genePositions(activeGenes, trnaAnnotation)
}

inactivePositions <- function (condition) {
    activeGeneNames <- getExpressedtRNAs(countsForCondition(condition))
    inactiveGeneNames <- setdiff(rownames(trnaUnfilteredRawCounts), activeGeneNames)
    inactiveGenes <- trnaUnfilteredRawCounts[inactiveGeneNames, ]
    genePositions(inactiveGenes, trnaUnfilteredAnnotation)
}

# Check whether active tRNA genes are more prone to have a certain mark close by
# than inactive tRNA genes
testActiveGeneMarkEnrichment <- function (active, inactive, marks, windowSize) {
    counts <- lapply(list(active, inactive),
                     countEnhancersPerGene, marks, windowSize)
    binarized <- lapply(counts, sign)
    binaryTable <- function (x)
        matrix(c(sum(x == 0), sum(x == 1)), ncol = 2)
    fisher.test(do.call(rbind, lapply(binarized, binaryTable)))
}

testMark <- function (condition, marks, windowSize)
    testActiveGeneMarkEnrichment(activePositions(condition),
                                 inactivePositions(condition),
                                 marks, windowSize)$p.value

if (! interactive()) {
    cat('# Generating colocalisation with H3K43me from Shen et al\n')

    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()

    loadMarks <- function (condition, mark) {
        x <- read.table(list.files(shenPath, sprintf('^%s.%s', condition, mark),
                                   full.names = TRUE),
                        col.names = c('Chr', 'Pos'), stringsAsFactors = FALSE)
        transform(x, Chr = sub('chr', '', Chr))
    }

    shenPath = '../common/data/shen-et-al/'
    marks <- c('h3k4me3', 'enhancer', 'h3k27ac', 'polII', 'ctcf')
    liverMarkConditions <- c('e14.5-liver', 'liver')

    liverMarks <- lapply(liverMarkConditions,
                         .(cond = setNames(lapply(marks, lp(loadMarks, cond)),
                                           marks))) %|%
        p(setNames, liverMarkConditions)

    testAllMarks <- function (windowSize)
        list(embryo = sapply(liverMarks$`e14.5-liver`,
                             .(mark = testMark('liver-e15.5', mark, windowSize))),
             adult = sapply(liverMarks$liver,
                            .(mark = testMark('liver-P29', mark, windowSize)))) %|%
        lp(map, unlist) %|% lp(do.call, rbind)

    windowSizes <- c(100, 500, 1000)
    allTests <- setNames(lapply(windowSizes, testAllMarks), windowSizes)

    mkdir('results/shen')

    writeExpressedData <- function (name, data) {
        filename = sprintf('results/shen/expressed-enrichment-pval-ws-%s.tsv',
                           name)
        write.table(data, filename, sep = '\t', quote = FALSE, col.names = NA)
    }

    map(writeExpressedData, names(allTests), allTests) %|% invisible
}
