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

#createRandomBackground <- function (acc, tissue) {
#    size <- nrow(subset(trnaAnnotation, Acceptor == acc))
#    trnas <- sample(rownames(trnaAnnotation), size)
#    cols <- tissueCols(tissue)
#    cor(t(trnaNormDataCond[trnas, cols]), method = 'spearman')
#}

createRotatedBackground <- function (acceptor, tissue) {
    geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    cols <- tissueCols(tissue)
    rotatedCor <- function (i)
        trnaNormDataCond[geneIds, …]
    lapply(indices(stages)[-1] - 1, rotatedCor)
}

createPermutedBackground <- function (acceptor, tissue) {
    geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    cols <- tissueCols(tissue)
    nbootstrap <- 100
    lapply(1 : nbootstrap,
           function (i) trnaNormDataCond[geneIds, cols[sample.int(length(cols))]])
}

createBackground <- function (acc, tissue)
    createPermutedBackground(acc, tissue)

plotBackgroundCorrelation <- function (acceptor)
    plotCorrelation(trnaBackgroundCor[[acceptor]],
                    main = paste('Acceptor background', acceptor))

regressClusters <- function(data, cluster1, cluster2)
    cor(colMeans(data[cluster2, ]), colMeans(data[cluster1, ]),
        method = 'spearman')

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

    allObs <- list()
    allBg <- list()
    ps <- list()
    clusterSizes <- list()

    findGeneClusters <- function (data) {
        if (nrow(data) < 2)
            return()

        pvclusters <- local({
            on.exit(sink())
            # Silence progress report.
            sink(file = '/dev/null')
            # pvclust can run into various errors. Nothing to be done.
            try(suppressWarnings(pvclust(t(data), nboot = 100)))
        })

        if (is(pvclusters, 'try-error'))
            return()
        clust <- pvpick(pvclusters)$clusters
        clusterSizes[[length(clusterSizes) + 1]] <- length(clust)
        clust
    }

    for (acceptor in unique(trnaAnnotation$Acceptor)) {
        cat(acceptor, '\n')
        geneIds <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
        cols <- tissueCols(tissue)
        data <- trnaNormDataCond[geneIds, cols]
        clust <- findGeneClusters(data)
        # Only consider cases with two clearly distinct clusters.
        if (length(clust) != 2) {
            warning('Skipping ', acceptor, ' because it has ',
                    length(clust), ' clusters', call. = FALSE)
            next
        }

        clusters <- do.call(rbind, lapply(clust, function (c) colMeans(data[c, ])))
        obs <- cor(clusters[1, ], clusters[2, ], method = 'spearman')
        permutations <- uperm(ncol(clusters))
        bg <- apply(permutations, ROWS,
                     function (p) cor(clusters[1, ], clusters[2, p],
                                      method = 'spearman'))
        p <- count(bg <= obs) / upermn(1 : ncol(clusters))
        ps[[acceptor]] <- p

        #obs <- regressClusters(data, clust[[1]], clust[[2]])
        #bg  <- sapply(trnaBackground, regressClusters, clust[[1]], clust[[2]])
        #allObs[[acceptor]] <- obs
        #allBg[[acceptor]] <- mean(bg)
    }
    ps <- unlist(ps)
    stripchart(ps, method = 'stack', pch = 20, xlab = 'p­values', xlim = c(0, 1),
               main = 'p­values for significant negative correlation\nof acceptor clusters')
    pOfH0 <- ks.test(ps, runif(1000))
    #' @TODO Add QQ plot instead of the silly stripchart

    haveSomeEvidence <- names(ps)[ps < 0.2]
    correlations <- unlist(sapply(haveSomeEvidence, function(n) trnaAcceptorCor[[n]][upper.tri(trnaAcceptorCor[[n]])]))
    hist(correlations, breaks = 25, freq = FALSE, col = 'grey', border = 'grey')
    lines(density(correlations), lwd = 2, col = colors[1])
}

# Visualise all correlation coefficients.
if (FALSE) {
    add <- FALSE
    for (acc in trnaAcceptorCor) {
        if (is.null(acc)) next
        hist(acc[upper.tri(acc)], breaks = 20, add = add, col = '#00000040',
             xlim = c(-1, 1), ylim = c(0, 12),
             main = 'Histogram of all correlation coefficients for all acceptors',
             xlab = 'Rank correlation coefficient of upper triangle of correlation matrix')
        add <- TRUE
    }
}
