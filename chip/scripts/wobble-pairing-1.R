reverseMap <- function (x)
    sapply(unique(unlist(codonsFor)), .(value = names(x)[grep(value, x)]))

# Possible codon-anticodon pairings

NONE <- character(0)

# Table copied manually from dos Reis & al., NAR (2004), Fig 1.
# This table is a nightmare; I have no idea how many typos remain – after hours
# of checking and removing numerous errors.
codonsFor <- list(
    AAA = NONE,
    GAA = c('TTT', 'TTC'),
    TAA = c('TTA', 'TTG'),
    CAA = 'TTG',

    AGA = c('TCT', 'TCC', 'TCA'),
    GGA = c('TCT', 'TCC'),
    TGA = c('TCA', 'TCG'),
    CGA = 'TCG',

    ATA = NONE,
    GTA = c('TAT', 'TAC'),
    TTA = NONE,
    CTA = NONE,

    ACA = NONE,
    GCA = c('TGT', 'TGC'),
    TCA = NONE,
    CCA = 'TGG',
# -----
    AAG = c('CTT', 'CTC', 'CTA'),
    GAG = c('CTT', 'CTC'),
    TAG = c('CTA', 'CTG'),
    CAG = 'CTG',

    AGG = c('CCT', 'CCC', 'CCA'),
    GGG = c('CCT', 'CCC'),
    TGG = c('CCA', 'CCG'),
    CGG = 'CCG',

    ATG = NONE,
    GTG = c('CAT', 'CAC'),
    TTG = c('CAA', 'CAG'),
    CTG = 'CAG',

    ACG = c('CGT', 'CGC', 'CGA'),
    GCG = c('CGT', 'CGC'),
    TCG = c('CGA', 'CGG'),
    CCG = 'CGG',
# -----
    AAT = c('ATT', 'ATC', 'ATA'),
    GAT = c('ATT', 'ATC', 'ATA'),
    TAT = 'ATA',
    CAT = 'ATG',

    AGT = c('ACT', 'ACC', 'ACA'),
    GGT = c('ACT', 'ACC'),
    TGT = c('ACA', 'ACG'),
    CGT = 'ACG',

    ATT = NONE,
    GTT = c('AAT', 'AAC'),
    TTT = c('AAA', 'AAG'),
    CTT = 'AAG',

    ACT = NONE,
    GCT = c('AGT', 'AGC'),
    TCT = c('AGA', 'AGG'),
    CCT = 'AGG',
# -----
    AAC = c('GTT', 'GTC', 'GTA'),
    GAC = c('GTT', 'GTC'),
    TAC = c('GTA', 'GTG'),
    CAC = 'GTG',

    AGC = c('GCT', 'GCC', 'GCA'),
    GGC = c('GCT', 'GCC'),
    TGC = c('GCA', 'GCG'),
    CGC = 'GCG',

    ATC = NONE,
    GTC = c('GAT', 'GAC'),
    TTC = c('GAA', 'GAG'),
    CTC = 'GAG',

    ACC = c('GGT', 'GGC', 'GGA'),
    GCC = c('GGT', 'GGC'),
    TCC = c('GGA', 'GGG'),
    CCC = 'GGG'
)

anticodonsFor <- reverseMap(codonsFor)

anticodonNames = unique(names(codonsFor))

fillAnticodonRows <- function (rows)
    let(mat = setNames(rows[anticodonNames], anticodonNames),
        `[<-`(mat, is.na(mat), 0))

wobbleWeights <- local({
    getWeights <- function (x) setNames(rep(1 / length(x), length(x)), x)
    withWeights <- lapply(anticodonsFor, getWeights)
    do.call(rbind, lapply(withWeights, fillAnticodonRows))
})

weightedWobblePairing <- function (codonWeights)
    colSums(wobbleWeights[names(codonWeights), ] * codonWeights)

mkdir('plots/wobble/all-alt')

local({
    mrna <- relativeData(codonUsageData)

    # The Python scripts counts stop codon usage. Remove them.
    mrna <- mrna[grep('TAA|TAG|TGA', rownames(mrna), invert = TRUE), ]

    mrna <- apply(mrna, COLS, weightedWobblePairing)

    trnaCodons <- groupby(trnaNormDataCond, trnaAnnotation[rownames(trnaNormDataCond), 'Acceptor'])
    # Unlike for the normal correlation, we do not reverse complement the
    # anticodons, since we are translating the codons instead.

    # But of course the *potential* anticodons contain anticodons which are not
    # utilised in mouse. We need to add these to the `trna` list.
    onlym <- setdiff(rownames(mrna), rownames(trnaCodons))
    trnaCodons[onlym, ] <- 0
    # And ensure that the row order is the same.
    trna <- relativeData(trnaCodons[rownames(mrna), ])

    on.exit(dev.off())
    pdf('plots/wobble/all-alt/correlations-all-genes.pdf', width = 7, height = 10,
        family = plotFamily)

    plotCodonAnticodonCorrelations(trna, mrna,
        title = 'Correlation accounting for\nwobble rules in %s %s')
})
