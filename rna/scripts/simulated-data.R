source('scripts/usage.R')

if (! interactive()) {
    cat('# Generating shuffled expression data tables\n')
    mrnaLoadData()
    mrnaNormalizeData()
    generateCodonUsageData()
    loadAminoAcids()
    loadCodonMap()
    resampleCodonUsage()

    prepare <- function (data) {
        data <- groupby(data, geneticCode[rownames(data), 'AA'])
        # Enforce uniform oder between tRNA and mRNA plots.
        data <- data[aminoAcids$Short, ]
        rownames(data) <- aminoAcids[aminoAcids$Short == rownames(data), 'Long']
        relativeData(data)
    }

    aaSampleMatrix <- map(prepare, codonSampleMatrix)
    aaSampleMeans <- do.call(cbind, map(rowMeans, aaSampleMatrix))

    aaMatrix <- relativeData(aaUsageData[, grep('liver', colnames(aaUsageData))])
    rownames(aaMatrix) <- aminoAcids[match(rownames(aaMatrix), aminoAcids$Short), 'Long']
    aaMatrix <- aaMatrix[aminoAcids$Long, ]

    mkdir('results/usage-sampling')

    map(.(data, name =
          write.table(data,
                      sprintf('results/usage-sampling/simulated-amino-acids-%s.tsv', name),
                      sep = '\t', quote = FALSE, col.names = FALSE)),
        aaSampleMatrix, names(aaSampleMatrix)) %|% invisible

    write.table(aaSampleMeans,
                'results/usage-sampling/simulated-amino-acids-means.tsv',
                sep = '\t', quote = FALSE, col.names = NA)

    write.table(aaMatrix, 'results/usage-sampling/real-amino-acids-means.tsv',
                sep = '\t', quote = FALSE, col.names = NA)
}
