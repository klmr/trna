source('../common/scripts/basic.R')
source('../common/scripts/perm.R')
source('scripts/load-data.R')
source('scripts/pvclust-patched.R')

library(gplots)

tissueCols <- function (tissue) {
    cols <- grep(tissue, colnames(trnaNormDataCond))
    cols[sapply(stages, grep, colnames(trnaNormDataCond)[cols])]
}

plotExpressionChange <- function (acceptor, tissue) {
    cols <- tissueCols(tissue)
    geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    if (length(geneIds) < 2)
        return()
    heatmap(as.matrix(trnaNormDataCond[geneIds, cols]), Colv = NA,
            col = progressColors,
            main = paste('Acceptor', acceptor, 'for', readable(tissue)),
            labCol = readable(colnames(trnaNormDataCond[, cols])))
}

isoacceptorCorrelations <- function (acceptor, tissue) {
    geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    cols <- tissueCols(tissue)
    if (length(geneIds) >= 2)
        cor(t(trnaNormDataCond[geneIds, cols]), method = 'spearman')
}

plotCorrelation <- function (data, ...)
    heatmap.2(data, symm = TRUE, col = contrastColors, Colv = TRUE,
              dendrogram = 'row', tracecol = last(colors), trace = 'none', ...)

plotAcceptorCorrelation <- function (acceptor)
    plotCorrelation(trnaAcceptorCor[[acceptor]], main = paste('Acceptor', acceptor))

scrambledCor <- function (x, method = 'pearson',
                          permutation = sample.int(nrow(x))) {
    result <- matrix(NA, nrow = ncol(x), ncol = ncol(x))

    for (i in  1 : (ncol(x) - 1))
        for (j in (i + 1) : ncol(x))
            result[i, j] <- cor(as.numeric(x[, i]),
                                as.numeric(x[permutation, j]),
                                method = method)

    result[lower.tri(result)] <- t(result)[lower.tri(result)]
    colnames(result) <- colnames(x)
    rownames(result) <- colnames(x)
    result
}

acceptorData <- function (acceptor, tissue)
    let(geneIds = rownames(subset(trnaAnnotation, Acceptor == acceptor)),
        cols = tissueCols(tissue),
        trnaNormDataCond[geneIds, cols])

findGeneClusters <- function (data) {
    if (nrow(data) < 2)
        return()

    # pvclust can run into various errors. Nothing to be done.
    pvclusters <- silent(try(suppressWarnings(pvclust(t(data), nboot = 100)),
                         silent = TRUE))

    if (! is(pvclusters, 'try-error'))
        pvpick(pvclusters)$clusters
}

clusterMeans <- function (acceptor, clust, data) {
    # Only consider cases with two clearly distinct clusters.
    if (length(clust) != 2) {
        message('Skipping ', acceptor, ' because it has ',
                length(clust), ' clusters')
        return()
    }
    do.call(rbind, lapply(clust, function (c) colMeans(data[c, ])))
}

compensationAnalysis <- function (tissue) {
    #' @TODO Is rank correlation really appropriate here?
    corMethod <- 'spearman'
    acceptors <- unique(trnaAnnotation$Acceptor)
    data <- setNames(map(acceptorData, acceptors, tissue), acceptors)
    clust <- map(findGeneClusters, data)
    clusters <- map(clusterMeans, names(clust), clust, data)
    clusters <- filter(neg(is.null), clusters)
    message('Using ', length(clusters), ' out of ', length(clust), ' isoacceptor families')

    observations <- map(.(c = cor(c[1, ], c[2, ], method = corMethod)), clusters)
    observations <- unlist(observations)
    permutations <- uperm(ncol(clusters[[1]]))
    background <- map(.(c = apply(permutations, ROWS,
                                  .(p = cor(c[1, ], c[2, p],
                                            method = corMethod)))), clusters)
    totalBackground <- setNames(do.call(c, background), NULL)
    test <- ks.test(observations, totalBackground, alternative = 'greater')

    structure(list(tissue = tissue,
                   background = totalBackground,
                   observations = observations,
                   test = test,
                   clusterSizes = table(unlist(map(length, clust)))),
              class = 'compensation')
}

print.compensation <- function (x, ...) {
    cat(sprintf('Frequencies of cluster sizes for analysis of %s:', x$tissue))
    print(x$clusterSizes)
    cat(sprintf('\np-value = %0.5f (ks.test)\n', x$test$p.value))
    x
}

plot.compensation <- function (x, ...) {
    hist(x$background, breaks = 25, col = 'grey', border = 'grey',
         main = 'Background distribution of correlations',
         xlab = 'Correlation coefficient',
         ylab = 'Frequency of correlation coefficient')
    invisible(map(.(x = abline(v = x, col = colors[1])), x$observations))
    par(usr = c(0, 1, 0, 1))
    text(1, 0.9, 'Observed\ncorrelations', col = colors[1], pos = 2)
    text(1, 0.1, bquote(italic(p) == .(sprintf('%0.3f', x$test$p.value))),
         col = colors[1], pos = 2)
}

plotHistAndDensity <- function (correlations, breaks = 25, ...) {
    hist(correlations, breaks = breaks, freq = FALSE, col = 'grey', border = 'grey',
         xlab = 'Correlation coefficient', ylim = c(0, 1))
    lines(density(correlations), lwd = 2, col = colors[1])
}

hist.compensation <- function (x, ...) {
    trnaAcceptorCor <- sapply(as.character(unique(trnaAnnotation$Acceptor)),
                              isoacceptorCorrelations, x$tissue)
    negCorAcceptors <- names(which(x$observations < -0.5))
    correlations <- map(.(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        negCorAcceptors) %|% unlist

    plotHistAndDensity(correlations, ...)
}

antihist <- function (x) {
    allAcceptors <-as.character(unique(trnaAnnotation$Acceptor))
    trnaAcceptorCor <- sapply(allAcceptors, isoacceptorCorrelations, x$tissue)
    negCorAcceptors <- names(which(x$observations < -0.5))
    use <- setdiff(allAcceptors, negCorAcceptors)
    use <- filter(.(x = ! is.null(trnaAcceptorCor[[x]])), use)
    correlations <- map(.(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        use) %|% unlist

    plotHistAndDensity(correlations)
}

chisq.test.compensation2 <- function (x) {
    bg <- table(x$background)
    obs <- table(x$observations)

    contingency <- rbind(as.vector(bg), rep(0, length(bg)))
    contingency[2, match(names(obs), names(bg))] <- as.vector(obs)
    chisq.test(contingency[2, ], p = contingency[1, ] / sum(contingency[1, ]))
}

chisq.test.default <- stats::chisq.test

chisq.test <- function (x, ...) UseMethod('chisq.test')

plot.compensation2 <- function (x, col = c('gray', 'red'), ...) {
    plot(density(x$background, adjust = 1.2), col = col[1], ...)
    lines(density(x$observations, adjust = 1.2), col = col[2], ...)
    text(0, 0.2, x$acceptor)
    test <- chisq.test(x)
    text(0, 0.1, bquote(italic(p) == .(sprintf('%0.2f', test$p.value))))
}

if (! interactive()) {
    cat('# Generating compensation analysis plots\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    mkdir('plots/compensation')

    # Plot two examples for publication.
    examples <- c('CAG', 'GCC')

    for (codon in examples) {
        map(.(tissue = {
            on.exit(dev.off())
            pdf(sprintf('plots/compensation/%s-%s.pdf', tissue, codon))
            plotExpressionChange(codon, tissue)
        }), tissues) %|% invisible
    }

    # Alternative analysis: pairwise gene correlations, no clusters.

    tissue <- 'liver'
    testSet <- unique(trnaAnnotation$Acceptor)

    compensationFile <- '../common/cache/compensation.RData'

    loaded <- tryCatch(load(compensationFile), error = .(. = ''))
    if (! identical(loaded, 'compensationData')) {
        message('Re-generating compensation analysis data')
        progress(0, length(testSet))
        compensationData <- setNames(map(.(acceptor, i = {
            correlations <- isoacceptorCorrelations(acceptor, tissue)
            if (is.null(correlations))
                return()
            correlations <- correlations[upper.tri(correlations)]

            geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
            cols <- tissueCols(tissue)
            genes <- trnaNormDataCond[geneIds, cols]

            permutations <- permn(ncol(genes))

            permCor <- function (genes, gene, perm)
                map(.(other = cor(as.numeric(genes[gene, ]),
                                  as.numeric(genes[other, perm]),
                                  method = 'spearman')),
                    (1 : nrow(genes))[-gene])

            map(.(gene = map(.(p = permCor(genes, gene, p)), permutations)),
                1 : nrow(genes)) -> background

            progress(i, length(testSet))
            structure(list(acceptor = acceptor,
                           observations = correlations,
                           background = unlist(background)),
                      class = 'compensation2')
        }), testSet, indices(testSet)), testSet)
        progress(1, 1)

        compensationData <- filter(neg(is.null), compensationData)

        save(compensationData, file = compensationFile)
    }

    local({
        on.exit(dev.off())
        pdf('plots/compensation/all-density.pdf', height = 10, width = 12, family = plotFamily)
        par(mfrow = c(8, 6), oma = rep(0, 4), mar = rep(0, 4))
        invisible(map(p(plot, col = c('gray', colors[1]), lwd = 2,
                        xlim = c(-1, 1), ylim = c(0, 1), bty = 'n', main = '',
                        xaxt = 'n', yaxt = 'n'),
                      compensationData))
    })

    map(.(codon = {
            on.exit(dev.off())
            pdf(sprintf('plots/compensation/evidence-%s.pdf', codon))
            plot(compensationData[[codon]], col = c('gray', colors[1]), lwd = 2,
                 xlim = c(-1, 1), ylim = c(0, 1),
                 xlab = 'Correlation coefficient',
                 main = paste('Compensation test for', codon))
        }), examples)

    pvalues <- unlist(map(item('p.value') %.% chisq.test, compensationData))
}
