source('../common/scripts/basic.R')
source('scripts/load-data.R')
source('scripts/de.R')

local({
    oldwd <- getwd(); on.exit(setwd(oldwd))
    setwd('../rna')
    source('scripts/de.R')
})

setdiff.data.frame <- function (x, y)
    x[setdiff(rownames(x), rownames(y)), ]

setdiff.default <- setdiff

setdiff <- function (x, y)
    UseMethod('setdiff')

colocalizationForContrast <- local({
    # Holds the already calculated differential expression test results for all
    # parameter combinations.
    cache <- list()

    # NOTE Of course Hadley did it first … package ‘memoise’.
    readCache <- function (data, tissue, a, b, threshold) {
        type <- as.character(substitute(data))
        key <- paste(type, tissue, a, b, threshold)
        result <- cache[[key]]
        if (is.null(result)) {
            result <- nbinomTest(data, paste(tissue, a, sep = '-'), paste(tissue, b, sep = '-'))
            cache[[key]] <- result
        }
        result
    }

    #' @param tissue the tissue
    #' @param a condition 1
    #' @param b condition 2
    #' @param threshold the threshold (FDR) at which to call significance
    #' @param windowSize The size of the window under consideration for neighboring genes
    function (tissue, a, b, threshold, windowSize) {
        warning('threshold unused')
        require(plyr) # for adply

        getDe <- function (data) subset(data, ! is.na(padj) & padj <= threshold)

        annotate <- function (data, annotation)
            annotation[data$id, c('Chr', 'Start', 'End')]

        prepare <- function (data, annotation) {
            all <- data[[tissue]][[a]][[b]]
            de <- getDe(all)
            up <- subset(de, foldChange < 1)
            down <- subset(de, foldChange > 1)
            no <- setdiff(all, de)
            lapply(list(all = all, de = de, up = up, down = down, no = no),
                   annotate, annotation)
        }

        countInNeighborhood <- function (trna, mrna)
            count(mrna$Chr == trna$Chr &
                  mrna$End >= trna$Start - windowSize &
                  mrna$Start <= trna$End + windowSize)

        if (! all(mrnaAnnotation$Start <= mrnaAnnotation$End))
            stop('Start <= End required for mrnaAnnotation!')

        trnaDe <- readCache(trnaCds, tissue, a, b, threshold)
        mrnaDe <- readCache(mrnaCds, tissue, a, b, threshold)
        trna <- prepare(trnaDe, trnaAnnotation)
        mrna <- prepare(mrnaDe, mrnaAnnotation)

        nearDe <- list(de = c(adply(trna$up, ROWS, countInNeighborhood, mrna$up)$V1,
                              adply(trna$down, ROWS, countInNeighborhood, mrna$down)$V1),
                       bg = c(adply(trna$up, ROWS, countInNeighborhood, mrna$all)$V1,
                              adply(trna$down, ROWS, countInNeighborhood, mrna$all)$V1))
        #    nearDe <- list(de = adply(trna$de, ROWS, countInNeighborhood, mrna$de)$V1,
        #                   bg = adply(trna$de, ROWS, countInNeighborhood, mrna$all)$V1)
        nearNonDe <- list(de = adply(trna$no, ROWS, countInNeighborhood, mrna$de)$V1,
                          bg = adply(trna$no, ROWS, countInNeighborhood, mrna$all)$V1)

        nearDe$ratio <- nearDe$de / nearDe$bg
        nearNonDe$ratio <- nearNonDe$de / nearNonDe$bg

        structure(list(de = nearDe$ratio,
                       non = nearNonDe$ratio,
                       tissue = tissue, a = a, b = b,
                       test = ks.test(nearDe$ratio, nearNonDe$ratio),
                       threshold = threshold, windowSize = windowSize),
                  class = 'colocalization')
    }
})

print.colocalization <- function (x, ...) {
    cat(sprintf('Colocalization between %s %s and %s\n\n',
                  readable(x$tissue), readable(x$a), readable(x$b)))
    cat(sprintf('D = %0.4f, p-value = %0.5f (ks.test)\n',
                test$statistic['D'], test$p.value))
    cat(sprintf('Parameters: ϑ = %s, w = %s\n', x$threshold, x$windowSize))
    x
}

plot.colocalization <- function (x, ...) {
    require(Hmisc)
    Ecdf(c(x$de, x$non),
         group = c(rep('DE', length(x$de)),
                   rep('non-DE', length(x$non))),
         main = sprintf('%s %s–%s (ϑ = %s, w = %s)',
                        readable(x$tissue), readable(x$a), readable(x$b),
                        x$threshold, x$windowSize),
         xlab = 'x', ylab = 'Proportion ≤ x', ...)
    #' @TODO Use proper math typesetting
    text(0.8, 0.1, sprintf('p = %0.2g', x$test$p.value))
}

ks.test.colocalization <- function (x)
    x$test

ks.test.default <- ks.test

ks.test <- function (x, y)
    UseMethod('ks.test')

testColocalization <- function () {
}

if (! interactive()) {
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    mrnaLoadData()
    mrnaSetupCountDataSet()
    mrnaPairwiseDifferentialExpression()
}

#' @TODO Everything below this is old, remove.

trnaLoadData()
trnaSetupCountDataSet()
trnaNormalizeData()
trnaDeAnalysis()

# Simple caching of results to speed up loading
mrnaContrasts <- new.env(parent = emptyenv())

# FIXME Ugly hack. Consolidate with other mRNA data.

mrnaLoadContrast <- function (tissue, stageA, stageB) {
    key <- paste(tissue, stageA, stageB)
    if (exists(key, mrnaContrasts))
        return(mrnaContrasts[[key]])

    stageA <- gsub('\\.', '', stageA)
    stageB <- gsub('\\.', '', stageB)
    dataFile <- file.path('../rna/data/deseq',
                          sprintf('%s_%s_VS_%s_%s.genes_de.tsv',
                                  tissue, stageA,
                                  tissue, stageB))
    mrnaContrasts[[key]] <- read.table(dataFile, header = TRUE, row.names = 1)
}

mrnaDeGenesForContrast <- function (contrast, threshold = 0.05)
    deGenesForContrast(contrast, threshold)

mrnaNonDeGenesForContrast <- function (contrast, threshold = 0.05)
    nonDeGenesForContrast(contrast, threshold)

#' `exclude` has the format list(list(Chr, From, To), list(Chr, From, To), ...)
testColocalization <- function (tissue, stageA, stageB, pval, windowSize, exclude) {
    stageBIndex <- grep(stageB, names(stageContrasts[[tissue]][[stageA]]))
    trnaContrast <- stageContrasts[[tissue]][[stageA]][[stageBIndex]]
    # Only use actively expressed genes.
    trnaContrast <- trnaContrast[trnaContrast$baseMeanA > 10 &
                                 trnaContrast$baseMeanB > 10, ]

    trnaGenes <- list(de = deGenesForContrast(trnaContrast, pval),
                      no = nonDeGenesForContrast(trnaContrast, pval))

    mrnaContrast <- mrnaLoadContrast(tissue, stageA, stageB)
    mrnaContrast$GOterm <- NULL
    mrnaContrast$GO <- NULL
    mrnaContrast$KEGG <- NULL
    mrnaContrast$lname <- NULL
    mrnaContrast <- mrnaContrast[mrnaContrast$source == 'protein_coding', ]
    mrnaContrast$source <- NULL
    mrnaContrast$chr <- regmatches(mrnaContrast$locus,
                                   regexpr('^[^:]+', mrnaContrast$locus))
    noChr <- regmatches(mrnaContrast$locus, regexpr('(\\d+)\\.{2}(\\d+)', mrnaContrast$locus))
    mrnaContrast$start <- as.numeric(regmatches(noChr, regexpr('^\\d+', noChr)))
    mrnaContrast$end <- as.numeric(regmatches(noChr, regexpr('\\d+$', noChr)))
    # Fix inverted start and end, we don't care about strands
    for (i in 1 : nrow(mrnaContrast)) {
        if (mrnaContrast[i, 'start'] > mrnaContrast[i, 'end']) {
            tmp <- mrnaContrast[i, 'start']
            mrnaContrast[i, 'start'] <- mrnaContrast[i, 'end']
            mrnaContrast[i, 'end'] <- tmp
        }
    }

    mrnaGenes <- list(de = mrnaDeGenesForContrast(mrnaContrast, pval),
                      no = mrnaNonDeGenesForContrast(mrnaContrast, pval))

    filterOutRegion <- function (genes, region) {
        filtered <- with(genes, Chr == region$Chr & Start >= region$Start & End <= region$End)
        cat(sprintf('Excluded %d genes\n', length(which(filtered))))
        genes[! filtered, ]
    }

    countInVicinity <- function (set, ann, width)
        nrow(subset(set, chr == ann$Chr &
                    end >= ann$Start - width &
                    end <= ann$End + width))

    counts <- list(de = matrix(nrow = nrow(trnaGenes$de), ncol = 2),
                   no = matrix(nrow = nrow(trnaGenes$no), ncol = 2))

    for (m in c('de', 'no')) {
        ann <- trnaAnnotation[trnaGenes[[m]]$id, ]
        ann$Chr <- sub('chr', '', regmatches(rownames(ann), regexpr('chr[^.]+', rownames(ann))))
        for (region in exclude)
            ann <- filterOutRegion(ann, region)

        for (i in 1 : nrow(trnaGenes[[m]])) {
            # Count number of DE mRNAs and total mRNAs in vicinity
            counts[[m]][i, 1] <- countInVicinity(mrnaGenes$de, ann[i, ], windowSize)
            counts[[m]][i, 2] <- countInVicinity(mrnaGenes$no, ann[i, ], windowSize) +
                                 counts[[m]][i, 1]
        }

        counts[[m]] <- cbind(counts[[m]], counts[[m]][, 1] / counts[[m]][, 2])
        colnames(counts[[m]]) <- c('DE', 'All', 'Ratio')
    }

    ks <- ks.test(counts$de[, 'Ratio'], counts$no[, 'Ratio'])

    require(Hmisc)
    Ecdf(c(counts$de[, 'Ratio'], counts$no[, 'Ratio']),
         main = sprintf('Distribution of ratio of proximal mRNAs to tRNAs\nin %s (%s vs %s) with p-val cutoff %s in %skb window',
                        tissue, stageA, stageB, pval, windowSize / 1000),
         xlab = 'Ratio (# diff. expr. mRNAs) / (# all mRNAs)', lwd = 2,
         group = c(rep(1, nrow(counts$de)), rep(2, nrow(counts$no))), col = c('black', 'blue'))
    text(0.8, 0.1, sprintf('p = %0.2f', ks$p.value))
    legend('topleft', legend = c('Diff. expr. tRNAs', 'Non-diff. expr. tRNAs'),
           fill = c('black', 'blue'), border = FALSE, bty = 'n')

    list(counts = counts, ks = ks)
}

plotColocalization <- function () {
    stageA <- 'e15.5'
    stageB <- 'P22'
    #exclude <- list(list(Chr = 13, Start = 21000000, End = 24000000))
    exclude <- list()

    for (tissue in c('liver', 'brain')) {
        for (pval in c(0.1, 0.05, 0.01)) {
            for (ws in c(10, 50, 100)) {
                path <- paste('plots/colocalization',
                              if (length(exclude) > 0) 'filtered', '', sep = '/')
                pdf(sprintf('%s/%s-%s-%s-pval_%s-window_%skb.pdf',
                            path, tissue, stageA, stageB, pval, ws), family = plotFamily)
                cat(sprintf('%s-%s-%s-pval_%s-window_%skb\n',
                            tissue, stageA, stageB, pval, ws))
                print(testColocalization(tissue, stageA, stageB, pval, ws * 1000, exclude))
                dev.off()
            }
        }
    }
}

mkdir('plots/colocalization/filtered')
plotColocalization()

if (FALSE) {
# Plot location of DE tRNAs {{{

# Data

chrSizes <- read.table('../common/data/mm9.sizes.tsv', row.names = 1)
chrSizesNames <- rownames(chrSizes)
chrSizes <- chrSizes[, 1]
names(chrSizes) <- chrSizesNames

trnaLoc <- trnaAnnotation[, c('Start', 'End')]
trnaLoc$Chr <- sub('\\..*', '', rownames(trnaAnnotation))

mrnaDeData <- list()

for (tissue in tissues) {
    for (stage in stages[-1]) {
        filename <- sprintf('../rna/data/deseq/%s_e155_VS_%s_%s.genes_de.tsv',
                            tissue, tissue, sub('\\.', '', stage))
        mrnaDeData[[tissue]][[paste('e15.5 vs', stage)]] <-
            read.table(filename, header = TRUE)
    }
}

# If we don't use `stringsAsFactors=F` then the following operation kills the RAM.
mrnaLocations <- read.table('../common/data/gene.annot.tsv',
                            stringsAsFactors = FALSE,
                            header = TRUE)[, c('ID', 'locus')]
mrnaChr <- lapply(mrnaLocations$locus, unlist %.% partial(strsplit, ':'))
mrnaLoci <- lapply(mrnaChr, unlist %.% partial(strsplit, '\\.\\.') %.% partial(slice, 2))

# This SEGFAULTs even though it works fine on smaller data!
#mrnaLocations$chr <- unlist(lapply(mrnaChr, lpartial(paste, 'chr', sep = '') %.% partial(slice, 1)))
mrnaChr <- unlist(lapply(mrnaChr, partial(slice, 1)))
mrnaLocations$Chr <- vapply(mrnaChr, function (x) paste('chr', x, sep = ''), '')
mrnaLocations$Start <- as.numeric(unlist(lapply(mrnaLoci, partial(slice, 1))))
mrnaLocations$End <- as.numeric(unlist(lapply(mrnaLoci, partial(slice, 2))))
rm(mrnaChr, mrnaLoci)
rownames(mrnaLocations) <- mrnaLocations$ID
mrnaLocations$ID <- NULL
mrnaLocations$locus <- NULL
#colnames(mrnaLocations) <- tolower(colnames(mrnaLocations))

# Analysis

getActive <- function (tissue)
    lapply(tissue, function (stage) stage[stage$baseMean > 20, 'id'])

allActivetRNAs <- c(getActive(allStages$liver), getActive(allStages$brain))
allActivetRNAs <- Reduce(union, allActivetRNAs)
active <- trnaLoc[allActivetRNAs, ]

isNum <- function (str)
    suppressWarnings(! is.na(as.numeric(str)))

toNum <- function (str)
    vapply(str, function (x) if (isNum(x)) as.numeric(x) else NA, 0)

chromosomes <- as.vector(names(chrSizes)[1:21])
chromosomes <- chromosomes[order(toNum(sub('chr', '', chromosomes)))]

plotChromosomeRegion <- function (chr, background, de, from = 1, to = chrSizes[chr], axis = FALSE) {
    background <- background[background$Chr == chr, ]
    de <- de[de$Chr == chr, ]
    graphics::plot(background$Start,
         rep(1, nrow(background)),
         xlim = c(from, to), ylim = c(0, 2), xaxs = 'i', axes = axis,
         xlab = 'Position', ylab = '',
         col = '#00000033', pch = 20)
    if (! axis)
        abline(0, 0)
    points(de$Start, rep(1.3, nrow(de)), col = '#FF000033', pch = 20)
    text(from, 1.8, chr, pos = 4)
}

detrnaGenesLoc <- function (tissue, stage) {
    stageDeData <- allStages[[tissue]][[paste('e15.5 vs', stage)]]
    trnaLoc[stageDeData[stageDeData$padj < 0.05 & ! is.na(stageDeData$padj), 'id'], ]
}

deGenesLoc <- function (tissue, stage) {
    stageDeData <- mrnaDeData[[tissue]][[paste('e15.5 vs', stage)]]
    mrnaLocations[stageDeData[stageDeData$padj < 0.05 & ! is.na(stageDeData$padj), 'id'], ]
}

for (tissue in tissues) {
    for (stage in stages[-1]) {
        de <- detrnaGenesLoc(tissue, stage)
        par(mfrow = c(length(chromosomes), 1), mar = c(0, 0, 0, 0))
        for (chr in chromosomes)
            plotChromosomeRegion(chr, active, de)
    }
}

inRange <- function (x, from, to)
    x[x > from & x < to]

par(mfrow = c(length(chromosomes), 1), mar = c(0, 0, 0, 0))
for (chr in chromosomes)
    plotChromosomeRegion(chr, active, detrnaGenesLoc('liver', 'P22'))

deMrna <- deGenesLoc('liver', 'P22')

# Region of interest
roi <- c(21000000, 24000000)
plotChromosomeRegion('chr13', active, detrnaGenesLoc('liver', 'P22'), roi[1], roi[2], axis = TRUE)
background <- mrnaLocations[mrnaLocations$Chr == 'chr13', 'start']
points(background, rep(0.5, length(background)), col = '#00000033', pch = 20)
points(deMrna$Start, rep(0.7, nrow(deMrna)), col = '#0000FF33', pch = 20)
dd <- density(inRange(deMrna$Start, roi[1], roi[2]), adjust = 0.1)
lines(dd$x, dd$y * 300000 + 0.7, col = 'blue')

library(zoo)

discreteDensity <- function (data, range, width, by = width / 2)
    rollapply(seq(range[1], range[2]), width = width, by = by,
              function (w) length(inRange(data, w[1], last(w))))

trnaDeLiverP22 <- detrnaGenesLoc('liver', 'P22')
trnaDeLiverP22chr13 <- trnaDeLiverP22[trnaDeLiverP22$chr == 'chr13', 'start']
numDeTrna <- discreteDensity(trnaDeLiverP22chr13, roi, 100000)
numAllTrna <- discreteDensity(active[active$chr == 'chr13', 'start'], roi, 100000)
ratioDe <- numDeTrna / numAllTrna
numDeMrna <- discreteDensity(deMrna$start, roi, 100000)

# Visualise the relative expressions accumulated in sliding windows:
{
    baseDens <- density(numAllTrna)
    plot(baseDens$x, baseDens$y, type = 'n')
    lines(baseDens, lwd = 2)
    lines(density(numDeTrna), col = 'red', lwd = 2)
    lines(density(numDeMrna), col = 'blue', lwd = 2)
}

# }}}
}
