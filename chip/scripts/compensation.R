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

plotBackgroundCorrelation <- function (acceptor)
    plotCorrelation(trnaBackgroundCor[[acceptor]],
                    main = paste('Acceptor background', acceptor))

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

if (! interactive()) {
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()

    #' @TODO Use as starting point to refactor in to function
    tissue <- 'liver'
    trnaAcceptorCor <- sapply(as.character(unique(trnaAnnotation$Acceptor)),
                              isoacceptorCorrelations, tissue)

    trnaBackgroundCor <- sapply(as.character(unique(trnaAnnotation$Acceptor)),
                                createRandomBackground, tissue)

    acceptorData <- function (acceptor)
        let(geneIds = rownames(subset(trnaAnnotation, Acceptor == acceptor)),
            cols = tissueCols(tissue),
            trnaNormDataCond[geneIds, cols])

    findGeneClusters <- function (data) {
        if (nrow(data) < 2)
            return()

        # pvclust can run into various errors. Nothing to be done.
        pvclusters <- silent(try(suppressWarnings(pvclust(t(data), nboot = 100))))

        if (! is(pvclusters, 'try-error'))
            pvpick(pvclusters)$clusters
    }

    clusterMeans <- function (clust, data) {
        # Only consider cases with two clearly distinct clusters.
        if (length(clust) != 2) {
            warning('Skipping ', acceptor, ' because it has ',
                    length(clust), ' clusters', call. = FALSE)
            return()
        }
        do.call(rbind, lapply(clust, function (c) colMeans(data[c, ])))
    }

    #' @TODO Is rank correlation really appropriate here?
    corMethod <- 'spearman'

    set.seed(123)

    acceptors <- unique(trnaAnnotation$Acceptor)
    data <- setNames(map(acceptorData, acceptors), acceptors)
    clust <- map(findGeneClusters, data)
    clusters <- map(clusterMeans, clust, data)
    clusters <- filter(neg(is.null), clusters)
    observations <- map(fun(c = cor(c[1, ], c[2, ], method = corMethod)), clusters)
    observations <- unlist(observations)
    background <- map(fun(c = apply(permutations, ROWS,
                                    function(p) cor(c[1, ], c[2, p],
                                                    method = corMethod))),
                      clusters)
    totalBackground <- do.call(c, background)
    test <- ks.test(observations, totalBackground)

    hist(totalBackground, breaks = 25, col = 'grey', border = 'grey',
         main = 'Background distribution of correlations',
         xlab = 'Correlation coefficient',
         ylab = 'Frequency of correlation coefficient')
    map(fun(x = abline(v = x, col = colors[1])), observations)
    par(usr = c(0, 1, 0, 1))
    text(1, 0.9, 'Observed\ncorrelations', col = colors[1], pos = 2)
    text(1, 0.1, bquote(italic(p) == .(sprintf('%0.3f', test$p.value))),
         col = colors[1], pos = 2)

    negCorAcceptors <- names(which(observations < -0.5))
    correlations <- map(fun(n = trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]),
                        negCorAcceptors) %|% unlist

    hist(correlations, breaks = 25, freq = FALSE, col = 'grey', border = 'grey')
    lines(density(correlations), lwd = 2, col = colors[1])
}
