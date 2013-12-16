source('../common/scripts/basic.R')
source('../common/scripts/perm.R')
source('scripts/load-data.R')
source('scripts/pvclust-patched.R')

library(gplots)

tissueCols <- function (tissue) {
    cols <- grep(tissue, colnames(trnaNormDataCond))
    cols[sapply(stages, grep, colnames(trnaNormDataCond)[cols])]
}

plotProgression <- function (acceptor, tissue) {
    cols <- tissueCols(tissue)
    geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    if (length(geneIds) < 2)
        return()
    heatmap.2(as.matrix(trnaNormDataCond[geneIds, cols]), Colv = FALSE,
              dendrogram = 'row', main = paste('Acceptor', acceptor),
              col = contrastColors, tracecol = last(colors), trace = 'none')
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

plotTestVis <- function (tissue, totalBackground, observations, test) {
    on.exit(dev.off())
    pdf(sprintf('plots/compensation/%s.pdf', tissue))
    hist(totalBackground, breaks = 25, col = 'grey', border = 'grey',
         main = 'Background distribution of correlations',
         xlab = 'Correlation coefficient',
         ylab = 'Frequency of correlation coefficient')
    invisible(map(fun(x = abline(v = x, col = colors[1])), observations))
    par(usr = c(0, 1, 0, 1))
    text(1, 0.9, 'Observed\ncorrelations', col = colors[1], pos = 2)
    text(1, 0.1, bquote(italic(p) == .(sprintf('%0.3f', test$p.value))),
         col = colors[1], pos = 2)
}

plotCorrelationDistribution <- function (tissue, trnaAcceptorCor, observations) {
    on.exit(dev.off())
    pdf(sprintf('plots/compensation/%s-correlations.pdf', tissue))
    negCorAcceptors <- names(which(observations < -0.5))
    correlations <- map(fun(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        negCorAcceptors) %|% unlist

    hist(correlations, breaks = 25, freq = FALSE, col = 'grey', border = 'grey')
    lines(density(correlations), lwd = 2, col = colors[1])
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
    totalBackground <- do.call(c, background)
    test <- ks.test(observations, totalBackground, alternative = 'greater')

    plotTestVis(tissue, totalBackground, observations, test)

    trnaAcceptorCor <- sapply(as.character(unique(trnaAnnotation$Acceptor)),
                              isoacceptorCorrelations, tissue)

    plotCorrelationDistribution(tissue, trnaAcceptorCor, observations)

    structure(list(clusterSizes = table(unlist(map(length, clust))),
                   p.value = test$p.value),
              class = 'compensation')
}

print.compensation <- function (x, ...) {
    cat('Frequencies of cluster sizes for analysis:')
    print(x$clusterSizes)
    cat(sprintf('\np-value = %0.5f (ks.test)\n', x$p.value))
    x
}

if (! interactive()) {
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    mkdir('plots/compensation')
    set.seed(123)
    map(compensationAnalysis, tissues)
}
