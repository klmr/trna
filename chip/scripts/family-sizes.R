# Test correlation between (expressed) isoacceptor family size and codon usage.

source('scripts/wobble-pairing-2.R')

rawIsoacceptorFamilySizes <- as.list(by(trnaUnfilteredAnnotation,
                                        trnaUnfilteredAnnotation$Acceptor,
                                        nrow)) %|% unlist

# Add nonexistent anticodons by selecting against preferred anticodons

isoacceptorFamilySizes <- rawIsoacceptorFamilySizes[preferredAnticodons]
trna <- do.call(cbind, map(.(. = isoacceptorFamilySizes), colnames(mrna)))
trna[duplicateAnticodons] <- trna[duplicateAnticodons] * codonUsageProp

plotIsoacceptorFamilySize <- function () {
    on.exit(dev.off())
    pdf('plots/wobble/only-missing-alt/correlations-family-sizes.pdf',
        width = 7, height = 10, family = plotFamily)
    plotCodonAnticodonCorrelations(trna, mrna,
                                   xlab = 'Anticodon isoacceptor family size')
}

# And, just for reference, do it without wobble pairing as well:

plotIsoacceptorFamilySizeNoWobbling <- function () {
    pairedAnticodons <- revcomp(rownames(codonUsageData))
    familySizeMatrix <- do.call(cbind, map(.(. = rawIsoacceptorFamilySizes[pairedAnticodons]),
                                           colnames(codonUsageData)))

    familySizeMatrix[is.na(familySizeMatrix)] <- 0

    on.exit(dev.off())
    pdf('plots/usage/codons-isoacceptor-family-size-scatter.pdf',
        width = 7, height = 10, family = plotFamily)
    plotCodonAnticodonCorrelations(familySizeMatrix, relativeData(codonUsageData),
                                   xlab = 'Anticodon isoacceptor family size',
                                   excludeZeros = TRUE)
}
