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

allDuplicated <- function (x)
    Reduce(`|`, lapply(x[duplicated(x)], p(`==`, x)))

# Calculate the codon usage proportions of all codon pairs competing for the
# same anticodon.

unadjustedCodonUsage <- lapply(alternativeAnticodons,
                               .(anticodon = mrna[preferredAnticodons == anticodon, ]))

codonUsageProp <- lapply(unadjustedCodonUsage, apply, COLS, .(x = x / sum(x)))
codonUsageProp <- do.call(rbind, codonUsageProp)

# Multiply anticodon usage of dual-use anticodons by codon usage proportions.

duplicateAnticodons <- allDuplicated(preferredAnticodons)

trna[duplicateAnticodons, ] <- trna[duplicateAnticodons, ] * codonUsageProp

plotCodonAnticodonCorrelations <- function (trna, mrna, excludeZeros = FALSE) {
    maxima <- findPlotMaxima(trna, mrna)
    par(mfrow = c(4, 3))

    map(.(tissue = {
        map(.(stage = {
            data <- columnsForCondition(trna, mrna, tissue, stage)
            colors <- if (excludeZeros)
                ifelse(data$trna == 0, last(colors), tissueColor[tissue]) else
                    tissueColor[tissue]
            plot(data$trna, data$mrna,
                 xlab = 'Proportion of tRNA isoacceptors',
                 ylab = 'Proportion of mRNA codon usage',
                 main = sprintf('Correlation accounting for\nwobble rules in %s %s', readable(stage), readable(tissue)),
                 xlim = c(0, maxima['trna']), ylim = c(0, maxima['mrna']),
                 col = colors, pch = 20, las = 1)
            if (excludeZeros)
                data <- data[data$trna != 0, ]
            model <- lm(mrna ~ trna, data)
            abline(model)
            par(usr = c(0, 1, 0, 1))
            rho <- cor(data$trna, data$mrna, method = 'spearman')
            r2 <- cor(data$trna, data$mrna, method = 'pearson')
            prho <- cor.test(data$trna, data$mrna, method = 'spearman')$p.value
            pr2 <- cor.test(data$trna, data$mrna, method = 'pearson')$p.value
            message(tissue, '-', stage, ': prho=', prho, ' pr2=', pr2)
            text(1, 0, bquote(atop(' ' ~ italic(p) == .(sprintf('%.2f', prho)) ~ (rho == .(sprintf('%.2f', rho))),
                                   italic(p) == .(sprintf('%.2f', pr2)) ~ (R^2 == .(sprintf('%.2f', r2))))),
                 adj = c(1.1, -0.1))

            rho
        }), stages) %|% unlist
    }), tissues)
}

mkdir('plots/wobble/only-missing-alt')

local({
    on.exit(dev.off())
    pdf('plots/wobble/only-missing-alt/correlations-all-genes.pdf', width = 7,
        height = 10, family = plotFamily)
    plotCodonAnticodonCorrelations(trna, mrna)
})
