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

# For DESeq3’s new result object
setdiff.DataFrame <- setdiff.data.frame

setdiff.default <- setdiff

setdiff <- function (x, y)
    UseMethod('setdiff')

#' @param tissue the tissue
#' @param a condition 1
#' @param b condition 2
#' @param threshold the threshold (FDR) at which to call significance
#' @param windowSize The size of the window under consideration for neighboring genes
colocalizationForContrast <- function (tissue, a, b, threshold, windowSize) {
    require(plyr) # for adply

    getDe <- function (data) subset(data, ! is.na(padj) & padj <= threshold)

    annotate <- function (data, annotation)
        let(selector = if (is.null(data$id)) rownames else item('id'),
            annotation[selector(data), c('Chr', 'Start', 'End')])

    prepare <- function (data, annotation) {
        all <- data[[tissue]][[a]][[b]]
        de <- getDe(all)
        # NB: Yes, `< 0` is correct: we look at contrast A–B and want to know which
        # genes are more highly expressed in contrast A, i.e. whose expression gets
        # lowered in B compared to A.
        up <- subset(de, log2FoldChange < 0)
        # Ditto
        down <- subset(de, log2FoldChange > 0)
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

    trna <- prepare(trnaDeResults, trnaAnnotation)
    mrna <- prepare(mrnaDeResults, mrnaAnnotation)

    nearDe <- list(de = adply(trna$up, ROWS, countInNeighborhood, mrna$up)$V1,
                   bg = adply(trna$up, ROWS, countInNeighborhood, mrna$all)$V1)
    nearNonDe <- list(de = adply(trna$no, ROWS, countInNeighborhood, mrna$up)$V1,
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

print.colocalization <- function (x, ...) {
    cat(sprintf('Colocalization between %s %s and %s\n\n',
                  readable(x$tissue), readable(x$a), readable(x$b)))
    cat(sprintf('D = %0.4f, p-value = %0.5f (ks.test)\n',
                x$test$statistic['D'], x$test$p.value))
    cat(sprintf('Parameters: ϑ = %s, w = %s, n = %s\n',
                x$threshold, x$windowSize, length(x$de)))
    x
}

plot.colocalization <- function (x, ...) {
    require(Hmisc)
    Ecdf(c(x$de, x$non),
         group = c(rep('DE', length(x$de)),
                   rep('non-DE', length(x$non))),
         main = substitute(tissue ~ a %<->% b ~
                          (paste(theta1 == t, ', ', italic(w) == ww)),
                           list(tissue = readable(x$tissue),
                                a = readable(x$a),
                                b = readable(x$b),
                                t = x$threshold,
                                ww = x$windowSize)),
         xlab = bquote(italic(x)), ylab = bquote('Proportion' <= italic(x)),
         label.curves = FALSE, subtitles = FALSE, ...)
    text(0.8, 0.1, bquote(italic(p) == .(x$test$p.value)))
    legend('topleft', legend = c('diff. expr.', 'non diff. expr.'), bty = 'n', ...)
}

ks.test.colocalization <- function (x)
    x$test

ks.test.default <- ks.test

ks.test <- function (x, y)
    UseMethod('ks.test')

testColocalization <- function () {
    # Test all parameter sets.
    stages <- list(liver = c('e15.5', 'P22'), brain = c('P4', 'P29'))
    path <- 'plots/colocalization'
    mkdir(path)
    maxProgress <- length(tissues) * 3 * 3
    currentProgress <- 0
    # Needed for progress counter
    env <- environment()

    pvalues <- lapply(tissues, function (tissue) {
        lapply(c(0.1, 0.05, 0.01), function (threshold) {
            lapply(c(10000, 50000, 100000), function (windowSize) {
                progress(currentProgress, maxProgress)
                assign('currentProgress', get('currentProgress', envir = env) + 1,
                       envir = env)

                stage <- stages[[tissue]]
                stat <- colocalizationForContrast(tissue, stage[1], stage[2],
                                                  threshold, windowSize)

                cond <- sprintf('%s-theta_%s-w_%s', tissue, threshold, windowSize)
                on.exit(dev.off())
                pdf(file.path(path, cond, ext = 'pdf'), family = plotFamily)
                plot(stat, col = colors, lwd = 3)
                ks.test(stat)$p.value
            })
        })
    })
    progress(currentProgress, maxProgress)
    pvalues
}

if (! interactive()) {
    cat('# Generating colocalization analysis plots\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaPairwiseDiffentialExpression()
    mrnaLoadData()
    mrnaPairwiseDifferentialExpression()

    pvalues <- testColocalization()
    adjusted <- unname(p.adjust(unlist(pvalues), method = 'bonferroni'))
    cat('Bonferroni-adjusted p-values of colocalization analyses:\n')
    map(p(cat, '\n'), adjusted) %|% invisible
}
