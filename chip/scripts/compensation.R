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

    observations <- map(fun(c = cor(c[1, ], c[2, ], method = corMethod)), clusters)
    observations <- unlist(observations)
    permutations <- uperm(ncol(clusters[[1]]))
    background <- map(fun(c = apply(permutations, ROWS,
                                    fun(p = cor(c[1, ], c[2, p],
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
    invisible(map(fun(x = abline(v = x, col = colors[1])), x$observations))
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
    correlations <- map(fun(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        negCorAcceptors) %|% unlist

    plotHistAndDensity(correlations, ...)
}

antihist <- function (x) {
    allAcceptors <-as.character(unique(trnaAnnotation$Acceptor))
    trnaAcceptorCor <- sapply(allAcceptors, isoacceptorCorrelations, x$tissue)
    negCorAcceptors <- names(which(x$observations < -0.5))
    use <- setdiff(allAcceptors, negCorAcceptors)
    use <- filter(fun(x = ! is.null(trnaAcceptorCor[[x]])), use)
    correlations <- map(fun(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        use) %|% unlist

    plotHistAndDensity(correlations)
}

if (! interactive()) {
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    mkdir('plots/compensation')
    set.seed(123)
    (compensationData <- map(compensationAnalysis, tissues))

    for (codon in c('CAG', 'AGT')) {
        map(fun(tissue = {
            on.exit(dev.off())
            pdf(sprintf('plots/compensation/%s-%s.pdf', tissue, codon))
            plotExpressionChange(codon, tissue)
        }), tissues) %|% invisible
    }

    map(fun(x = {
        on.exit(dev.off())
        pdf(sprintf('plots/compensation/%s-correlations.pdf', x$tissue),
            width = 8, height = 4, family = plotFamily)
        par(mfrow = c(1, 2))
        hist(x)
        antihist(x)
    }), compensationData) %|% invisible

    map(fun(x = {
        on.exit(dev.off())
        pdf(sprintf('plots/compensation/%s.pdf', tissue), family = plotFamily)
        plot(x)
    }), compensationData) %|% invisible

    codons <- setNames(map(fun(x = names(x$observations)[x$observations < 0]),
                           compensationData), NULL)
    differentCodons <- do.call(setdiff, codons)
    sharedCodons <- do.call(intersect, codons)
    table(geneticCode[differentCodons, ])
    table(geneticCode[sharedCodons, ])
}
