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

testUpregulatedGeneMarkEnrichment <- function (genes, marks, windowSize, threshold) {
    de <- deGenes(genes, threshold)
    # NB: Yes, `< 0` is correct: we look at contrast A–B and want to know which
    # genes are more highly expressed in contrast A, i.e. whose expression gets
    # lowered in B compared to A.
    up <- subset(de, log2FoldChange < 0)
    #down <- subset(de, log2FoldChange > 0)
    down <- setdiff(genes, up)

    testActiveGeneMarkEnrichment(genePositions(up, trnaAnnotation),
                                 genePositions(down, trnaAnnotation),
                                 marks, windowSize)$p.value
}

if (! interactive()) {
    cat('# Generate colocalisation with histone marks from Shen et al\n')

    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()

    cat('# Test for enrichment of marks near active tRNA genes\n')

    loadMarks <- function (condition, mark) {
        x <- read.table(list.files(shenPath, sprintf('^%s.%s', condition, mark),
                                   full.names = TRUE),
                        col.names = c('Chr', 'Pos'), stringsAsFactors = FALSE)
        transform(x, Chr = sub('chr', '', Chr))
    }

    shenPath = '../common/data/shen-et-al/'
    markNames <- c('h3k4me3', 'enhancer', 'h3k27ac', 'polII', 'ctcf')
    liverMarkConditions <- c('e14.5-liver', 'liver')
    brainMarkConditions <- c('e14.5-brain', 'cortex', 'cerebellum')

    loadAllMarks <- function (conditions)
        lapply(conditions,
               .(cond = lapply(markNames, lp(loadMarks, cond)) %|%
                 p(setNames, markNames))) %|%
        p(setNames, conditions)

    liverMarks <- loadAllMarks(liverMarkConditions)
    names(liverMarks) <- c('embryo', 'adult')
    brainMarks <- loadAllMarks(brainMarkConditions)
    brainMarks$embryo <- brainMarks$`e14.5-brain`
    brainMarks$adult <- map(rbind, brainMarks$cortex, brainMarks$cerebellum)

    marks = list(liver = liverMarks, brain = brainMarks)

    testAllExpressedMarks <- function (tissue, windowSize)
        list(embryo = sapply(marks[[tissue]]$embryo,
                             .(mark = testMark(sprintf('%s-e15.5', tissue),
                                               mark, windowSize))),
             adult = sapply(marks[[tissue]]$adult,
                            .(mark = testMark(sprintf('%s-P29', tissue),
                                              mark, windowSize)))) %|%
        lp(map, unlist) %|% lp(do.call, rbind)

    windowSizes <- c(100, 500, 1000)
    expressedTests <- map(.(tissue = map(lp(testAllExpressedMarks, tissue),
                                         windowSizes) %|% p(setNames, windowSizes)),
                          names(marks))

    mkdir('results/shen')

    writeExpressedData <- function (tissue, windowSize, data) {
        filename = sprintf('results/shen/pvalues-expressed-%s-ws-%s.tsv',
                           tissue, windowSize)
        write.table(data, filename, sep = '\t', quote = FALSE, col.names = NA)
    }

    map(.(tissue, data = map(.(ws, data = writeExpressedData(tissue, ws, data)),
                             names(data), data)),
        names(expressedTests), expressedTests) %|% invisible

    cat('# Test for enrichment of marks near differentially expressed tRNA genes\n')

    testAllUpregulatedMarks <- function (tissue, windowSize, threshold) {
        trnas <- trnaDeResults[[tissue]]$e15.5$P29
        test <- function (mark)
            testUpregulatedGeneMarkEnrichment(trnas, mark, windowSize, threshold)

        list(embryo = sapply(marks[[tissue]]$embryo, test),
             adult = sapply(marks[[tissue]]$adult, test)) %|%
        lp(map, unlist) %|% lp(do.call, rbind)
    }

    upregulatedTests <- map(.(tissue = map(lp(testAllUpregulatedMarks, tissue),
                                           windowSizes, 0.01) %|% p(setNames, windowSizes)),
                            names(marks))

    writeUpregulatedData <- function (tissue, windowSize, data) {
        filename = sprintf('results/shen/pvalues-de-%s-ws-%s.tsv',
                           tissue, windowSize)
        write.table(data, filename, sep = '\t', quote = FALSE, col.names = NA)
    }

    map(.(tissue, data = map(.(ws, data = writeUpregulatedData(tissue, ws, data)),
                             names(data), data)),
        names(upregulatedTests), upregulatedTests) %|% invisible
}
