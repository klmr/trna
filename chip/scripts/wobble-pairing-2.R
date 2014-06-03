# Find out which codons have no corresponding anticodon.
# We parsimoniously assume that only those codons require wobble pairing
# to form a codon/anticodon bond.

codonPairingPartners <- revcomp(as.character(unique(trnaAnnotation$Acceptor)))
# Remove stop codons.
allCodons <- setdiff(rownames(codonUsageData), c('TAA', 'TAG', 'TGA'))
unpairedCodons <- setdiff(allCodons, codonPairingPartners)

# Find alternative anticodons.
# The following code relies on the assumption that unpaired codons only have a
# single anticodon which recognises them through wobble pairing.
# This isn’t generally true, but it holds in mouse, and it makes things easier.

alternativeAnticodons <- anticodonsFor[unpairedCodons]
# And from that, subtract again those anticodons which do not exist in mice.
alternativeAnticodons <- lapply(alternativeAnticodons,
                                lp(filter, neg(p(`%in%`, unpairedCodons)
                                               %.% revcomp))) %|% unlist

preferredAnticodons <- revcomp(allCodons)
preferredAnticodons[unpairedCodons] <- alternativeAnticodons

mrna <- relativeData(codonUsageData)
# The Python scripts counts stop codon usage. Remove them.
mrna <- mrna[grep('TAA|TAG|TGA', rownames(mrna), invert = TRUE), ]

trna <- groupby(trnaNormDataCond, trnaAnnotation$Acceptor, sum)
trna <- relativeData(trna[preferredAnticodons, ])

# Calculate the codon usage proportions of all codon pairs competing for the
# same anticodon.

unadjustedCodonUsage <- lapply(alternativeAnticodons,
                               .(anticodon = mrna[preferredAnticodons == anticodon, ]))

codonUsageProp <- lapply(unadjustedCodonUsage, apply, COLS, .(x = x / sum(x)))
codonUsageProp <- do.call(rbind, codonUsageProp)

# Multiply anticodon usage of dual-use anticodons by codon usage proportions.

duplicateAnticodons <- allDuplicated(preferredAnticodons)

trna[duplicateAnticodons, ] <- trna[duplicateAnticodons, ] * codonUsageProp

mkdir('plots/wobble/only-missing-alt')

local({
    on.exit(dev.off())
    pdf('plots/wobble/only-missing-alt/correlations-all-genes.pdf', width = 7,
        height = 10, family = plotFamily)
    plotCodonAnticodonCorrelations(trna, mrna)
})
