source('scripts/correlation.R')

if (! interactive()) {
    cat('# Generating shuffled expression data tables\n')
    trnaLoadData()
    trnaSetupCountDataSet()
    trnaNormalizeData()
    trnaGroupFamilyAndType()
    loadAminoAcids()
    loadGeneticCode()
    resampleAcceptorAbundance()

    makeSameOrder <- p(`[`, aminoAcids$Long, quote(expr = ))

    isotypeSampleMatrix <- map(p(groupby, trnaUnfilteredAnnotation$Type) %|>%
                               makeSameOrder %|>%  relativeData,
                               acceptorSampleMatrix)
    isotypeSampleMeans <- do.call(cbind, map(rowMeans, isotypeSampleMatrix))

    onlyLiver <- trnaNormDataCond[, grep('liver', colnames(trnaNormDataCond))]
    isotypeMatrix <- onlyLiver[, vapply(stages, p(grep, colnames(onlyLiver)), numeric(1))] %|%
        p(groupby, trnaAnnotation$Type) %|% makeSameOrder %|% relativeData

    mkdir('results/usage-sampling')

    map(.(data, name =
          write.table(data,
                      sprintf('results/usage-sampling/simulated-isotypes-%s.tsv', name),
                      sep = '\t', quote = FALSE, col.names = FALSE)),
        isotypeSampleMatrix, names(isotypeSampleMatrix)) %|% invisible

    write.table(isotypeSampleMeans,
                'results/usage-sampling/simulated-isotypes-means.tsv',
                sep = '\t', quote = FALSE, col.names = NA)

    write.table(isotypeMatrix, 'results/usage-sampling/real-isotypes-means.tsv',
                sep = '\t', quote = FALSE, col.names = NA)
}
