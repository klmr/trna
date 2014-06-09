# Address reviewer's comment that co-regulated tRNA genes in anticodon iso-
# acceptor family, which show evidence of compensation with other genes in the
# same family, may cluster spatially on genome

source('scripts/de.R')

trnaLoadData()
trnaSetupCountDataSet()
trnaNormalizeData()
trnaPairwiseDiffentialExpression()

trnaClusterFile <- '../common/data/all_tRNA_genes_in_clusters.tsv'

trnaClusters <- read.table(trnaClusterFile, skip = 1, stringsAsFactors = FALSE,
                           sep = '\t', header = TRUE)
rm(trnaClusterFile)
rownames(trnaClusters) <- sub('_', '.', trnaClusters$ID)
trnaClusters <- with(trnaClusters,
                     data.frame(ID = clusters, Chr = chr.1,
                                Start = start.1, End = end.1,
                                Active = rownames(trnaClusters) %in% rownames(trnaAnnotation),
                                row.names = rownames(trnaClusters),
                                stringsAsFactors = FALSE))

genesInCluster <- by(trnaClusters, trnaClusters$ID, rownames)

areClustersUniformlyActive <- by(trnaClusters, trnaClusters$ID,
                                 function (x) all(x$Active) || all(! x$Active))

areClustersUniformlyActive <- setNames(as.logical(areClustersUniformlyActive),
                                       dimnames(areClustersUniformlyActive)[[1]])

nonUniformExpressedClusters <-
    names(areClustersUniformlyActive)[! areClustersUniformlyActive]

# Check whether either all tRNAs are differentially expressed in the same direction,
# or whether they are all non-DE.
# Ensure that at least one gene in the cluster is actually expressed, otherwise we
# include all-inactive clusters, making our number of coserved clusters artificially
# high.

liverDeGenes <- Reduce(union, lapply(trnaDeGenes, rownames))

trnaCheckDE <- function (trnas)
    let(genesAreDe = rownames(trnas) %in% liverDeGenes,
        atLeastOneActive = length(intersect(rownames(trnas),
                                            rownames(trnaAnnotation))) > 0,
        atLeastOneActive && (all(genesAreDe) || all(! genesAreDe)))

areClustersUniformlyDE <- by(trnaClusters, trnaClusters$ID, trnaCheckDE)
areClustersUniformlyDE <- setNames(as.logical(areClustersUniformlyDE),
                                   dimnames(areClustersUniformlyDE)[[1]])

nonUniformDeClusters <-
    names(areClustersUniformlyDE)[! areClustersUniformlyDE]
# These include clusters without any active gene. Filter these out.
hasActiveGenes <- .(x = nrow(subset(trnaClusters, ID == x & Active)) > 0)
nonUniformDeClusters <- filter(hasActiveGenes, nonUniformDeClusters)

trnasInNonUniformDeClusters <- rownames(subset(trnaClusters,
                                               ID %in% nonUniformDeClusters))

# ... and now what?

# Look at the isoacceptor compensation data to compare the isoacceptors that
# show compensation with those that don't show compensation.

# First, we compare the percentage of genes within each isoacceptor family that
# fall into clusters with a random background.

# For each gene, determine whether it is in a cluster with at least one other
# gene of the same isoacceptor family. Return the percentage of these genes.
isoacceptorGenesInClusterPercent <- function (acceptor, families)
    let(clusters = trnaClusters[families[[acceptor]], 'ID'],
        count(allDuplicated(clusters, NA)) / length(clusters))

isoacceptorGenes <- function (acc)
    rownames(subset(trnaAnnotation, Acceptor == acc))

allIsoacceptors <- unique(as.character(trnaAnnotation$Acceptor))
# Only isoacceptors for which we performed compensation analysis
compIsoacceptors <- filter(.(acc = nrow(subset(trnaAnnotation,
                                               Acceptor == acc)) > 5),
                           allIsoacceptors)

allFamilies <- sapply(allIsoacceptors, isoacceptorGenes)
compFamilies <- sapply(compIsoacceptors, isoacceptorGenes)

allPercentInCluster <- sapply(allIsoacceptors, isoacceptorGenesInClusterPercent,
                              allFamilies)
compPercentInCluster <- sapply(compIsoacceptors,
                               isoacceptorGenesInClusterPercent, compFamilies)

# Simulate random isoacceptor assignments and calculate numbers

familySizes <- sapply(allFamilies, length)

set.seed(536723) # Resulted from keyboard mashing

# Take N samples without replacement of sizes given by the second argument.
# Sizes must be nonzero.
splitSample <- function (x, sizes, prob = NULL) {
    splitChunks <- unlist(sapply(1 : length(sizes), .(i = rep(i, sizes[i]))))
    setNames(split(sample(x, sum(sizes), prob = prob), splitChunks),
             names(sizes))
}

randomBackground <- sapply(1 : 1000, .(i = {
    randomFamilies <- splitSample(rownames(trnaAnnotation), familySizes)
    percentInCluster <- sapply(names(randomFamilies),
                               isoacceptorGenesInClusterPercent, randomFamilies)
}))

local({
    on.exit(dev.off())
    pdf('plots/compensation/isoacceptor-enrichment-of-genomic-clusters.pdf',
        family = plotFamily)

    hist(randomBackground[, 1], xlim = c(0, 1), ylim = c(0, 25),
         col = transparent(last(colors), 0.15), border = NA, las = 1,
         xlab = 'Percentage of genes in clusters per isoacceptor',
         main = 'Predominance of clustered genes')

    apply(randomBackground[, -1], COLS,
          .(x = hist(x, col = transparent(last(colors)), border = NA,
                     add = TRUE)))

    hist(allPercentInCluster, col = transparent(colors[1]), border = NA,
         add = TRUE)

    legend('topright', bty = 'n', fill = c(last(colors), colors[1]), border = NA,
           legend = c('Background from random sampling', 'Isoacceptor families'))
})

# Test whether the (mean) genes-in-cluster percentage per isoacceptor family is
# higher than by chance (i.e. in the randomly generated data):

meanIsSmaller <- mean(allPercentInCluster) < colMeans(randomBackground)

# Test whether most of the means are smaller (H0: 50% smaller, 50% larger)
pClusterPercentageEqual <- binom.test(count(meanIsSmaller),
                                      length(meanIsSmaller))$p.value

message('p(H0): means are equal = ', pClusterPercentageEqual)

# Test whether isoacceptors with disproportionately high cluster coverage are
# enriched in set of isoacceptors showing evidence of compensation effects.

test <- names(allPercentInCluster[allPercentInCluster > mean(allPercentInCluster)])
compensatedIsoacceptors <- rownames(testValues[testValues$adjusted < 0.05, ])
uncompensatedIsoacceptors <- rownames(testValues[testValues$adjusted >= 0.05, ])

enrichedInCompensated <- count(test %in% compensatedIsoacceptors)
enrichedInUncompensated <- count(test %in% uncompensatedIsoacceptors)

message(enrichedInCompensated, ' out of ', length(compensatedIsoacceptors),
        sprintf(' (%2.0f%%)', 100 * enrichedInCompensated /
                length(compensatedIsoacceptors)),
        ' isoacceptors(compensated) are clustered more than average')
message(enrichedInUncompensated, ' out of ', length(uncompensatedIsoacceptors),
        sprintf(' (%2.0f%%)', 100 * enrichedInUncompensated /
                length(uncompensatedIsoacceptors)),
        ' isoacceptors(uncompensated) are clustered more than average')

# Next, we count the number of clusters used in an isoacceptor, to test the
# hypothesis that non-compensated isoacceptors use more clusters - in other
# words, that compensated isoacceptors' genes tend to cluster together.
# For this purposes, genes not in clusters count as distinct clusters.
# Cluster counts are normalised by isoacceptor family size, since bigger
# families have potentially more clusters.

clusterCount <- function (acceptor) {
    genes <- rownames(subset(trnaAnnotation, Acceptor == acceptor))
    clusters <- trnaClusters[genes, 'ID']
    count(! duplicated(clusters, incomparables = NA)) / length(genes)
}

compensatedClusters <- sapply(compensatedIsoacceptors, clusterCount)
uncompensatedClusters <- sapply(uncompensatedIsoacceptors, clusterCount)

percentages <- sort(union(compensatedClusters, uncompensatedClusters))
counttab <- sapply(percentages, .(p = c(count(compensatedClusters == p),
                                        count(uncompensatedClusters == p))))
significant <- chisq.test(counttab)
message('Difference between compensated and uncompensated isoacceptor ',
        'families: p = ', significant$p.value)

local({
    on.exit(dev.off())
    pdf('plots/compensation/cluster-enrichment-in-isoacceptors.pdf',
        family = plotFamily)
    hist(uncompensatedClusters, xlim = c(0, 1), col = colors[2],
         main = 'Distribution of counts of clusters in isoacceptor families',
         xlab = 'Count of clusters / size of family')
    hist(compensatedClusters, xlim = c(0, 1), col = transparent(colors[1]),
         add = TRUE)
    p <- wilcox.test(compensatedClusters, uncompensatedClusters)$p.value
    legend('topleft', legend = c('Clusters without compensation',
                                 'Clusters with compensation'),
           fill = c(colors[2], transparent(colors[1])), bty = 'n',
           border = NA)
    text(0.2, 1, sprintf('p = %0.3f', p))
})
